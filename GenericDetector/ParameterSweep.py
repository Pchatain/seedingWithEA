#!/usr/bin/env python3
from deap import cma
from deap import algorithms
from deap import tools
from deap import creator
from deap import base

from collections import OrderedDict

import pathlib
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import random
import subprocess
import multiprocessing
import numpy as np
import array

# User Definitions
PT_CUTS = [0.1, 1, 2.5, 3]
PT_CUT_COLORS = ["b", "r", "g", "k"]
assert len(PT_CUTS) > 0
assert len(PT_CUTS) <= len(PT_CUT_COLORS)
# Define bounds on parameters during training
MINS = OrderedDict([('maxPt', 0), ('impactMax', 0.1), ('deltaRMin', 0.25), ('sigmaScattering', 0.2), ('deltaRMax', 55), ('maxSeedsPerSpM', 0), ('radLengthPerSeed', 0.001)]) #, 0.001, 400, 5]
MAXS = OrderedDict([('maxPt', 12345), ('impactMax', 20), ('deltaRMin', 50), ('sigmaScattering', 10), ('deltaRMax', 150), ('maxSeedsPerSpM', 4), ('radLengthPerSeed', 0.02)])
# MAXS = [12345, 20, 50, 10, 150, 4, 0.02] #, 0.003, 600, 10]
# Dictionary of normalization coefficients
# because update for each parameter is drawn from the same normal distribution

NAME_TO_FACTOR = OrderedDict([('maxPt', 12000), ('impactMax', 1.0), ('deltaRMin', 5.0), ('sigmaScattering', 2.0), ('deltaRMax', 60.0), ('maxSeedsPerSpM', 1.0), ('radLengthPerSeed', 0.005)])
# myGuess = [12000, 20, 2.25, 3, 220, 0.95, 0.005] # ldmx
myGuess = [12000, 9, 1, 2.25, 60, 0.99, 0.005] # ttbar bad
NAME_TO_INDEX = {}
for i, oneName in enumerate(NAME_TO_FACTOR):
    # myGuess[i] *= 1.0 / NAME_TO_FACTOR[oneName]
    # MINS[oneName] *= 1.0 / NAME_TO_FACTOR[oneName]
    # MAXS[oneName] *= 1.0 / NAME_TO_FACTOR[oneName]
    NAME_TO_INDEX[oneName] = i
NSWEEPS = 5 # Number of times to sweep over all parameters
NEVALS = 30 # Number of different combinations to try 

plotDirectory = "ySweep_7params_scored1_ttbar_" + str(NSWEEPS) + "sweeps_" + str(NEVALS) + "evals_pop50_srange0.01-0.3_eval1" # Where to save the plots
plotDirectory += "/"
ttbarSampleInput = ['--input-dir', 'sim_generic_ATLASB_ttbar_e4_pu200_eta2.5/']
ttbarSampleBool = True


# Tags to use for reading output from seeding algorithm
mlTag = 'mlTag'
effTag, dupTag, seedsTag, truthTag = 'eff', 'dup', 'seeds', 'true' 


# Create the names with the parameters for the seeding algo
def createNamesAndParams(individual):
    names = []
    params = []
    for i, name in enumerate(NAME_TO_FACTOR):
        names.append(name)
        if (name == 'maxSeedsPerSpM'):
            params.append(int(individual[i]))
        else:
            params.append(individual[i])
    return names, params


# Print an individual
def indPrint(ind):
    newInd = normalizedToActual(ind)
    commandLineInd = ''
    for argName in newInd:
        commandLineInd += ("--sf-" + argName + " " +
                           str(newInd[argName]) + " ")
    print(commandLineInd)

# Format the input for the seeding algorithm. 
# Assumes program is in same directory as seeding algorithm
def paramsToInput(params, names):
    ret = ['./ActsExampleTestSeedAlgorithm',
           '--response-file', 'config_seeding_ml']
    if (ttbarSampleBool):
        ret.append(ttbarSampleInput[0])
        ret.append(ttbarSampleInput[1])
    if len(params) != len(names):
        raise Exception("Length of Params must equal names in paramsToInput")
    i = 0
    for param in params:
        arg = "--sf-" + names[i]
        ret.append(arg)
        paramValue = param
        ret.append(str(paramValue))
        i += 1
    return ret

# Opens a subprocess that runs the seeding algorithm and retrieves output using grep
# Returns efficiency, fake rate, and duplicate rate as percentages
def executeAlg(arg):
    p2 = subprocess.Popen(
        arg, bufsize=4096, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1_out, p1_err = p2.communicate()
    p1_out = p1_out.decode()
    p1_out = p1_out.strip().encode()
    p2 = subprocess.Popen(
        ['grep', mlTag], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output = p2.communicate(input=p1_out)[0].decode().strip()
    tokenizedOutput = output.split('\n')
    eventCount = {'dup': 0, 'eff': 0, 'seeds': 0, 'tSeeds': 0}
    scoreTotals = {'efficiency': -1, 'fakeRate': -1, 'duplicateRate': -1}
    # loop over all events
    for line in tokenizedOutput:
        if (line.find(mlTag) == -1):
            # print(f"word not found: {line}")
            continue
        # print(f"{line}")
        # Get data from one event
        seedingOutput = {'dup': -1, 'eff': -1, 'seeds': -1, 'tSeeds': -1}
        splittedLine = line.split(',')
        for word in splittedLine:
            if (word.find(dupTag) != -1):
                seedingOutput['dup'] = word[len(dupTag):]
                eventCount['dup'] += 1
            elif (word.find(effTag) != -1):
                seedingOutput['eff'] = (word[len(effTag):])
                eventCount['eff'] += 1
            elif (word.find(seedsTag) != -1):
                seedingOutput['seeds'] = (word[len(seedsTag):])
                eventCount['seeds'] += 1
            elif (word.find(truthTag) != -1):
                seedingOutput['tSeeds'] = (word[len(truthTag):])
                eventCount['tSeeds'] += 1
        if seedingOutput['eff'] == -1 or seedingOutput['seeds'] == -1 or seedingOutput['tSeeds'] == -1 or seedingOutput['dup'] == -1:
            continue
        efficiency, nSeeds, nTrueSeeds, nDup = float(seedingOutput['eff']), int(seedingOutput['seeds']), int(seedingOutput['tSeeds']), int(seedingOutput['dup'])
        if (nSeeds == 0):
            continue
        fakeRate = 100 * (nSeeds - nTrueSeeds) / nSeeds
        if (nTrueSeeds == 0):
          duplicateRate = 0
        else:
          duplicateRate = 100 * nDup / nTrueSeeds
        scoreTotals["efficiency"] += efficiency
        scoreTotals["fakeRate"] += fakeRate
        scoreTotals["duplicateRate"] += duplicateRate
    nEvents = eventCount['dup']
    # print(f"N events is {nEvents}")
    if nEvents != eventCount['eff'] or nEvents != eventCount['seeds'] or nEvents != eventCount['tSeeds']:
        print(f"Houston, we have a problem. One line was missing information. we had {nEvents} dup events")
        print("we passed in:")
        print(arg)
        print()
    avgScores = {}
    for score in scoreTotals.keys():
        if nEvents == 0:
            avgScores[score] = 0
        else:
            avgScores[score] = scoreTotals[score] / nEvents
    # print(avgScores)
    return avgScores

# Evaluates an individual and calculates a score
# First value - score - is used to rank individuals. Other 3 values considered only when first one is equal.
def evaluate(individual):
    names, params = createNamesAndParams(individual)
    arg = paramsToInput(params, names)
    r = executeAlg(arg)
    eff, fakeRate, duplicateRate = r['efficiency'], r['fakeRate'], r['duplicateRate']
    # MAX_SEEDS = 20000
    # seedsScore = 10 * float(seeds) / MAX_SEEDS
    # if (nTrueSeeds == 0):
    #     return -1.0, 100.0, 100.0
    effScore = (1 / (1 - (eff / 100)))
    penalty = fakeRate * duplicateRate / (1000) # min(effScore, 200) - penalty
    #effWeighted = eff
    effWeighted = eff - penalty
    return effWeighted, eff, fakeRate, duplicateRate


# Plots a graph of the best individual according to one of the 3 metrics and one parameter
# HOF = Hall Of Fame, records best indiviual over totallity of training so far
def plotHOF(scoreName, scores, paramName, paramValues):
    saveDir = plotDirectory + "hof/" + scoreName + "/"
    pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True)
    plotName = saveDir + paramName + ".png"
    fig, ax1 = plt.subplots()
    gen = range(1, NGEN + 1)
    line1 = ax1.plot(gen, scores, "b-", label="best " + scoreName)
    ax1.set_xlabel("Generation")
    ax1.set_ylabel(scoreName, color="b")
    for tl in ax1.get_yticklabels():
        tl.set_color("b")

    ax2 = ax1.twinx()
    line2 = ax2.plot(gen, paramValues, "r-", label=paramName)
    ax2.set_ylabel(paramName, color="r")
    for tl in ax2.get_yticklabels():
        tl.set_color("r")

    lns = line1 + line2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="best")
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close(fig)


# Creates a plot with the efficiency for different pT ranges
def plotPtRange(pTCutData):
    pathlib.Path(plotDirectory).mkdir(parents=True, exist_ok=True)
    plotName = plotDirectory + "Efficiency_vs_pT_ranges.png"
    fig, ax = plt.subplots()
    rangeToLine = {}
    for i, cutRange in enumerate(pTCutData):
        lineName = str(cutRange[0]) + "-" + str(cutRange[1]) + "GeV"
        rangeToLine[cutRange] = ax.plot(pTCutData[cutRange], PT_CUT_COLORS[i] + "-o", label=lineName)
    ax.axis([0, NGEN - 1, -1, 100])
    # lns = rangeToLine[pTCutData[(pT)]]
    # labs = [l.get_label() for l in lns]
    ax.legend(loc="best")
    plt.tight_layout()
    print("Now we should be saving the figure...")
    plt.savefig(plotName)
    plt.close(fig)


# Plots a graph of the best individual according to one of the 3 metrics and one parameter
# HOF = Hall Of Fame, records best indiviual over totallity of training so far
def plotEffFakeDup(sweepNum, effs, fakes, dups):
    saveDir = plotDirectory
    pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True)
    plotName = saveDir + "sweepNum_" + str(sweepNum) +  "effFakeDup.png"
    fig, ax1 = plt.subplots()
    gen = range(1, len(NAME_TO_FACTOR) + 1)
    line1 = ax1.plot(gen, effs, "b-", label="Efficiency")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Efficiency", color="b")
    for tl in ax1.get_yticklabels():
        tl.set_color("b")

    ax2 = ax1.twinx()
    line2 = ax2.plot(gen, fakes, "r-", label="Fake Rate")
    line3 = ax2.plot(gen, dups, "r--", label="Duplicate Rate")
    ax2.set_ylabel("Fake And Duplicate Rate", color="r")
    for tl in ax2.get_yticklabels():
        tl.set_color("r")

    lns = line1 + line2 + line3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="best")
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close(fig)

# Plot just the efficiency, fake rate and duplicate for the best individal (HOF) over time
def plotScores(efficiencies, fakeRates, duplicateRates):
    pathlib.Path(plotDirectory).mkdir(parents=True, exist_ok=True)
    plotName = plotDirectory + "Efficiency_FakeRate_DuplicateRate.png"
    # metrics = ["Efficiency", "FakeRate", "DuplicateRate"]
    fig, ax = plt.subplots()
    l1 = ax.plot(efficiencies, "b-", label="efficiency")
    l2 = ax.plot(fakeRates, "r-", label="fakeRate")
    l3 = ax.plot(duplicateRates, "g-", label="duplicateRate")
    ax.axis([0, NSWEEPS * len(NAME_TO_FACTOR) - 1, -1, 100])
    lns = l1 + l2 + l3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="best")
    plt.tight_layout()
    # print("Now we should be saving the figure..")
    plt.savefig(plotName)
    plt.close(fig)


creator.create("Fitness", base.Fitness, weights=(1.0, 1.0, -1.0, -1.0))
creator.create("Individual", array.array, typecode="d",
               fitness=creator.Fitness)

# def sortInds(individuals, fit_attr="fitness"):
#     return sorted(individuals, key=attrgetter(fit_attr), reverse=True)

def initPopulation(pcls, ind_init, myGuess):
    pop = pcls(ind_init(myGuess) for i in range(NEVALS))
    # for ind in pop:
    #     ind.strategy = scls(random.uniform(SMIN_INITIAL, SMAX_INITIAL) for _ in range(len(ind)))
    return pop

toolbox = base.Toolbox()
toolbox.register("evaluate", evaluate)
toolbox.register("sortAll", tools.selBest, k=NEVALS)
toolbox.register("population_guess", initPopulation, list,
                 creator.Individual)


# Create multiprocessing pool
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)


def createIndividuals(nEvals, mins, maxs, bestInd, oneName):
    pop = toolbox.population_guess(bestInd)
    # print(f"param is {oneName}")
    # print(f"max is {maxs[oneName]} min is {mins[oneName]}")
    paramValueRange = maxs[oneName] - mins[oneName]
    # print(f"param Value range is {paramValueRange}")
    stepSize = paramValueRange / nEvals
    # print(f"step size is {stepSize}")
    for i, ind in enumerate(pop):
        ind[NAME_TO_INDEX[oneName]] = mins[oneName] + i * stepSize
    return pop

def oneSweep(nEvals, mins, maxs, bestInd):
    effs, fakes, dups, scores = [], [], [], []
    for oneName in NAME_TO_FACTOR:
        pop = createIndividuals(nEvals, mins, maxs, bestInd, oneName)
        # print("pop is ")
        # print(pop)
        fitnesses = toolbox.map(toolbox.evaluate, pop)
        # print(fitnesses)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit
        # sortedPop = sorted(pop, key=attrgetter(fit_attr), reverse=True)
        newBestInd = toolbox.sortAll(pop)[0]
        iName = NAME_TO_INDEX[oneName]
        bestInd[iName] = newBestInd[iName]
        print(f"Best val for {oneName} is {newBestInd[iName]} with fitness {newBestInd.fitness.values}")
        score, eff, fake, dup = newBestInd.fitness.values
        scores.append(score)
        effs.append(eff)
        fakes.append(fake)
        dups.append(dup)
    bestConfig = []
    for val in bestInd:
        bestConfig.append(val)
    return bestConfig, scores, effs, fakes, dups

def main(nSweeps, nEvals, mins, maxs):
    # Initialize dictionary to keep track of Hall Of Fame (HOF) data
    hofData = {}
    for oneName in NAME_TO_FACTOR:
        hofData[oneName] = []
    # Keep track of HOF over time, so that we can evaluate efficiency over different pT ranges
    bestIndividuals = []
    # Data for plotting efficiency as a function of pT
    pTCutData = {}
    for i, cut in enumerate(PT_CUTS):
        cutRange = 0, 100
        if i == len(PT_CUTS) - 1:
            cutRange = cut, 12345
        else:
            cutRange = cut, PT_CUTS[i+1]
        pTCutData[cutRange] = []
    # Objects to keep track of seeding algorithm metrics for HOF
    scores = []
    efficiencies = []
    fakeRateList = []
    dupRateList = []
    bestInd = [] # best individual found after each sweep through parameters
    for val in myGuess:
        bestInd.append(val) # initialize best ind
    bestScore = 0
    for isweep in range(nSweeps):
        bestInd, scores, effs, fakes, dups = oneSweep(nEvals, mins, maxs, bestInd)
        newBestScore = max(scores)
        # print(f"Best ind is {bestInd} with {newBestScore} score")
        plotEffFakeDup(isweep, effs, fakes, dups)
        if (newBestScore - bestScore == 0):
          break
        else:
          bestScore = newBestScore
        
main(NSWEEPS, NEVALS, MINS, MAXS)
print(f"Wrote plots to {plotDirectory}")
