// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFilter.hpp"

#include <cmath>
#include <numeric>
#include <type_traits>

#include <fstream>
#include <iostream>
#include <sstream>

namespace Acts {

template <typename external_spacepoint_t, typename platform_t>
Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(std::move(config)) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);
}

template <typename external_spacepoint_t, typename platform_t>
template <typename sp_range_t>
std::vector<Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
          size_t print1 = 0;
    size_t print2 = 0;
  std::vector<Seed<external_spacepoint_t>> outputVec;
  for (auto spM : middleSPs) {
    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();

    // bottom space point
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>
        compatBottomSP;

    for (auto bottomSP : bottomSPs) {
      float rB = bottomSP->radius();
      float deltaR = rM - rB;
      // if r-distance is too big, try next SP in bin
      if (deltaR > m_config.deltaRMax) {
        continue;
      }
      // if r-distance is too small, continue because bins are NOT r-sorted
      if (deltaR < m_config.deltaRMin) {
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = (zM - bottomSP->z()) / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        continue;
      }
      compatBottomSP.push_back(bottomSP);
    }
    // no bottom SP found -> try next spM
    if (compatBottomSP.empty()) {
      continue;
    }

    std::vector<const InternalSpacePoint<external_spacepoint_t>*> compatTopSP;

    for (auto topSP : topSPs) {
      float rT = topSP->radius();
      float deltaR = rT - rM;
      // this condition is the opposite of the condition for bottom SP
      if (deltaR < m_config.deltaRMin) {
        continue;
      }
      if (deltaR > m_config.deltaRMax) {
        continue;
      }

      float cotTheta = (topSP->z() - zM) / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        continue;
      }
      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        continue;
      }
      compatTopSP.push_back(topSP);
    }
    if (compatTopSP.empty()) {
      continue;
    }
    // contains parameters required to calculate circle with linear equation
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;
    transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
    transformCoordinates(compatTopSP, *spM, false, linCircleTop);

    // create vectors here to avoid reallocation in each loop
    std::vector<const InternalSpacePoint<external_spacepoint_t>*> topSpVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSpM;
    size_t numBotSP = compatBottomSP.size();
    size_t numTopSP = compatTopSP.size();


    for (size_t b = 0; b < numBotSP; b++) {
      auto lb = linCircleBottom[b];
      float Zob = lb.Zo;
      float cotThetaB = lb.cotTheta;
      float Vb = lb.V;
      float Ub = lb.U;
      float ErB = lb.Er;
      float iDeltaRB = lb.iDeltaR;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1. + cotThetaB * cotThetaB);
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
      // scattering
      // but to avoid trig functions we approximate cot by scaling by
      // 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at
      // eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *=
          m_config.sigmaScattering * m_config.sigmaScattering;

      // clear all vectors used in each inner for loop
      topSpVec.clear();
      curvatures.clear();
      impactParameters.clear();
      for (size_t t = 0; t < numTopSP; t++) {
        auto lt = linCircleTop[t];

        // add errors of spB-spM and spM-spT pairs and add the correlation term
        // for errors on spM
        float error2 = lt.Er + ErB +
                       2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
                           iDeltaRB * lt.iDeltaR;

        float deltaCotTheta = cotThetaB - lt.cotTheta;
        float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
        float error;
        float dCotThetaMinusError2;
        // if the error is larger than the difference in theta, no need to
        // compare with scattering
        if (deltaCotTheta2 - error2 > 0) {
          deltaCotTheta = std::abs(deltaCotTheta);
          // if deltaTheta larger than the scattering for the lower pT cut, skip
          error = std::sqrt(error2);
          dCotThetaMinusError2 =
              deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
          // avoid taking root of scatteringInRegion
          // if left side of ">" is positive, both sides of unequality can be
          // squared
          // (scattering is always positive)

          if (dCotThetaMinusError2 > scatteringInRegion2) {
            continue;
          }
        }

        // protects against division by 0
        float dU = lt.U - Ub;
        if (dU == 0.) {
          continue;
        }
        // A and B are evaluated as a function of the circumference parameters
        // x_0 and y_0
        float A = (lt.V - Vb) / dU;
        float S2 = 1. + A * A;
        float B = Vb - A * Ub;
        float B2 = B * B;
        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        if (S2 < B2 * m_config.minHelixDiameter2) {
                    if (print1 < 1) {
          //std::cout << "Momentum cut " << 2 * m_config.pTPerHelixRadius / std::sqrt(B2/S2) << std::endl;
           print1++; }
          continue;
        } else if (print2 < 1) {
          //std::cout << "Momentum kept " << 2 * m_config.pTPerHelixRadius / std::sqrt(B2/S2) << std::endl; 
          print2++; 
        }
        // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
        float iHelixDiameter2 = B2 / S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2scatter = 4 * iHelixDiameter2 * m_config.pT2perRadius;
        // TODO: include upper pT limit for scatter calc
        // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
        // from rad to deltaCotTheta
        float p2scatter = pT2scatter * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if ((2 * m_config.pTPerHelixRadius * std::sqrt(S2) < m_config.maxPt * B) && (deltaCotTheta2 - error2 > 0) &&
            (dCotThetaMinusError2 >
             p2scatter * m_config.sigmaScattering * m_config.sigmaScattering)) {
          continue;
        } 
        // A and B allow calculation of impact params in U/V plane with linear
        // function
        // (in contrast to having to solve a quadratic function in x/y plane)
        float Im = std::abs((A - B * rM) * rM);

        if (Im <= m_config.impactMax) {
          topSpVec.push_back(compatTopSP[t]);
          // inverse diameter is signed depending if the curvature is
          // positive/negative in phi
          curvatures.push_back(B / std::sqrt(S2));
          impactParameters.push_back(Im);
        }
      }
      if (!topSpVec.empty()) {
        std::vector<std::pair<
            float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
            sameTrackSeeds;
        sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(
            *compatBottomSP[b], *spM, topSpVec, curvatures, impactParameters,
            Zob));
        seedsPerSpM.insert(seedsPerSpM.end(),
                           std::make_move_iterator(sameTrackSeeds.begin()),
                           std::make_move_iterator(sameTrackSeeds.end()));
      }
    }
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);
  }
  return outputVec;
}

template <typename external_spacepoint_t, typename platform_t>
void Seedfinder<external_spacepoint_t, platform_t>::transformCoordinates(
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec) const {
  float xM = spM.x();
  float yM = spM.y();
  float zM = spM.z();
  float rM = spM.radius();
  float varianceZM = spM.varianceZ();
  float varianceRM = spM.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  for (auto sp : vec) {
    float deltaX = sp->x() - xM;
    float deltaY = sp->y() - yM;
    float deltaZ = sp->z() - zM;
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY * sinPhiM;
    float y = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ * iDeltaR * bottomFactor;
    // VERY frequent (SP^3) access
    LinCircle l;
    l.cotTheta = cot_theta;
    // location on z-axis of this SP-duplet
    l.Zo = zM - rM * cot_theta;
    l.iDeltaR = iDeltaR;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    l.U = x * iDeltaR2;
    l.V = y * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    l.Er = ((varianceZM + sp->varianceZ()) +
            (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
           iDeltaR2;
    linCircleVec.push_back(l);
  }
}
}  // namespace Acts
