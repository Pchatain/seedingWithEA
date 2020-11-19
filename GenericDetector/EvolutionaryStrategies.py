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
import json
import array


# User Definitions
PT_CUTS = [0.1, 0.5, 1, 10]
PT_CUT_COLORS = ["b", "r", "g", "k"]
assert len(PT_CUTS) > 0
assert len(PT_CUTS) <= len(PT_CUT_COLORS)
# Define bounds on parameters during training
MINS = [1200, 0.1, 0.25, 0.2, 50, 0, 0.001] #, 0.001, 400, 5]
MAXS = [1234567, 20, 30, 10, 200, 4, 0.02] #, 0.003, 600, 10]
# Dictionary of normalization coefficients
# because update for each parameter is drawn from the same normal distribution
# NAME_TO_FACTOR = {'maxPt': 12000, 'impactMax': 1.0, 'deltaRMin': 5.0, 'sigmaScattering': 2.0, 'deltaRMax': 60.0, 'maxSeedsPerSpM': 1.0, 'radLengthPerSeed': 0.005}
NAME_TO_FACTOR = OrderedDict([('maxPt', 12000), ('impactMax', 1.0), ('deltaRMin', 5.0), ('sigmaScattering', 2.0), ('deltaRMax', 60.0), ('maxSeedsPerSpM', 1.0), ('radLengthPerSeed', 0.005)])
# NAME_TO_FACTOR = {'maxPt': 12000, 'impactMax': 1,
#                   'deltaRMin': 5, 'sigmaScattering': 2, 'deltaRMax': 60.0, 'collisionRegionMin': -3, 'collisionRegionMax': 3, 'maxSeedsPerSpM': 1, 'radLengthPerSeed': 0.005, 'bFieldInZ': 0.002, 'minPt': 500, 'cotThetaMax': 7.40627}
# myGuess = [12000, 1, 1, 2.25, 60, 0.95, 0.005] # good
myGuess = [25000, 5, 5, 0.2, 120, 0.95, 0.005] # bad
# myGuess = [12000, 10, 1, 2.25, 60, 0.99, 0.005] # ttbar
# myGuess = [1234, 9, 4.0, 2.0, 80, 0.9, 0.005]
NAME_TO_INDEX = {}
for i, oneName in enumerate(NAME_TO_FACTOR):
    myGuess[i] *= 1.0 / NAME_TO_FACTOR[oneName]
    MINS[i] *= 1.0 / NAME_TO_FACTOR[oneName]
    MAXS[i] *= 1.0 / NAME_TO_FACTOR[oneName]
    NAME_TO_INDEX[oneName] = i
NGEN = 2 # Number of generations
plotDirectory = "yGARBAGIO_________5" #"yES_7params_scored4_muon_gen" + str(NGEN) + "_pop50_srange0.01-0.3_eval1" # "zES_7params_scored1_muon_gen200_pop50_srange0.01-0.3_eval2" # Where to save the plots
plotDirectory += "/"
ttbarSampleInput = ['--input-dir', 'sim_generic_ATLASB_ttbar_e4_pu200_eta2.5/']
ttbarSampleBool = False

# Evolitionary Algorithm Parameters
NPOP = 50 # Population size
TournamentSize = 3 # Parameter used for selection

CXPB, MUTPB, SIGMA, INDPB = 0.5, 0.3, 0.1, 0.2
SMIN_INITIAL, SMAX_INITIAL = 0.3, 0.5 # max and minimum strategy used for the first mutation (separates the initial population which has all values the same)
smin, smax = 0.01, 0.1 # Minimum and maximum "strategies" to use for updating parameters. 
# a strategy specifies the standard deviation of the gaussin to draw mutations
# Each individual has a different strategy, and each parameter within an individual has a different strategy

# Tags to use for reading output from seeding algorithm
mlTag = 'mlTag'
effTag, dupTag, seedsTag, truthTag = 'eff', 'dup', 'seeds', 'true'

def normalizedToActual(normalizedInd):
    newDict = {}
    for i, name in enumerate(NAME_TO_FACTOR):
        if (name == 'maxSeedsPerSpM'):
            newDict[name] = int(NAME_TO_FACTOR[name] * normalizedInd[i])
        else:
            newDict[name] = NAME_TO_FACTOR[name] * normalizedInd[i]
    return newDict

# Create the names with the parameters for the seeding algo
def createNamesAndParams(individual):
    names = []
    params = []
    for i, name in enumerate(NAME_TO_FACTOR):
        names.append(name)
        if (name == 'maxSeedsPerSpM'):
            params.append(int(individual[i] * NAME_TO_FACTOR[name]))
        else:
            params.append(individual[i] * NAME_TO_FACTOR[name])
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
        arg, bufsize=4096, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
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
        fakeRate = 100 * (nSeeds - nTrueSeeds) / nSeeds
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
            print(f"no events found!")
            break
        avgScores[score] = scoreTotals[score] / nEvents
    # print(avgScores)
    return avgScores


creator.create("Fitness", base.Fitness, weights=(1.0, 1.0, -1.0, -1.0))
creator.create("Individual", array.array, typecode="d",
               fitness=creator.Fitness, strategy=None)
creator.create("Strategy", array.array, typecode="d")


def initPopulation(pcls, scls, ind_init, myGuess):
    pop = pcls(ind_init(myGuess) for i in range(NPOP))
    for ind in pop:
        ind.strategy = scls(random.uniform(SMIN_INITIAL, SMAX_INITIAL) for _ in range(len(ind)))
    return pop


toolbox = base.Toolbox()

# toolbox.register("attr_int", random.randint, INT_MIN, INT_MAX)
# toolbox.register("attr_flt", random.uniform, FLT_MIN, FLT_MAX)
# toolbox.register("individual", tools.initRepeat, creator.Individual,
#                  toolbox.attr_flt, n=N_CYCLES)
toolbox.register("population_guess", initPopulation, list, creator.Strategy,
                 creator.Individual, myGuess)

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

# Forces individual to stay within bounds after an update
def checkBounds(mins, maxs):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    if child[i] > maxs[i]:
                        child[i] = maxs[i]
                    elif child[i] < mins[i]:
                        child[i] = mins[i]
            return offspring
        return wrapper
    return decorator

# Makes sure strategy is within bounds
def checkStrategy(minstrategy, maxstrategy):
    def decorator(func):
        def wrappper(*args, **kargs):
            children = func(*args, **kargs)
            for child in children:
                for i, s in enumerate(child.strategy):
                    if s < minstrategy:
                        child.strategy[i] = minstrategy
                    elif s > maxstrategy:
                        child.strategy[i] = maxstrategy
            return children
        return wrappper
    return decorator


toolbox.register("evaluate", evaluate)
# toolbox.register("mate", tools.cxESBlend, alpha=0.5)
toolbox.register("mutate", tools.mutESLogNormal, c=1, indpb=INDPB)
# toolbox.register("mate", tools.cxTwoPoint) # We don't use crossover

# Mutate function draws from normal distribution with mean mu and std sigma
# Each parameter within an indibidual is mutated with indpb chance
# toolbox.register("mutate", tools.mutGaussian,
#                  mu=0.0, sigma=SIGMA, indpb=INDPB)  
# Population is selected by drawing tournsize individuals at random, and keeping the best one, NPOP times
toolbox.register("select", tools.selTournament, tournsize=TournamentSize)

# toolbox.decorate("mate",  checkBounds(MINS, MAXS))
toolbox.decorate("mutate", checkBounds(MINS, MAXS))
toolbox.decorate("mutate", checkStrategy(smin, smax))


def printStats(fits, length, title):
    mean = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    std = abs(sum2 / length - mean**2)**0.5

    print("  {0} Min {1}".format(title, min(fits)))
    print("  {0} Max {1}".format(title, max(fits)))
    print("  {0} Avg {1}".format(title, mean))
    print("  {0} Std {1}".format(title, std))

# Create multiprocessing pool
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)

# Plots a graph of the population according to one of the 3 metrics and one parameter
def plotLogbook(scoreName, paramName, logbook):
    saveDir = plotDirectory + "population/" + scoreName + "/"
    pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True) 
    plotName = saveDir + paramName + ".png"
    fig, ax1 = plt.subplots()
    gen = logbook.chapters[paramName].select("gen")
    # ax1.set_xticks(range(len(gen) + 1))
    scoresMax = logbook.chapters[scoreName].select("max")
    scoresMin = logbook.chapters[scoreName].select("min")
    line1 = ax1.plot(gen, scoresMax, "b-", label="max" + scoreName)
    line5 = ax1.plot(gen, scoresMin, "b--", label="min" + scoreName)
    ax1.set_xlabel("Generation")
    ax1.set_ylabel(scoreName, color="b")
    for tl in ax1.get_yticklabels():
        tl.set_color("b")

    ax2 = ax1.twinx()
    param_max = logbook.chapters[paramName].select("max")
    param_min = logbook.chapters[paramName].select("min")
    param_avg = logbook.chapters[paramName].select("avg")
    for i in range(len(param_max)):
        param_max[i] *= NAME_TO_FACTOR[paramName]
        param_min[i] *= NAME_TO_FACTOR[paramName]
        param_avg[i] *= NAME_TO_FACTOR[paramName]
    line2 = ax2.plot(gen, param_max, "r-", label="max")
    line3 = ax2.plot(gen, param_avg, "r:", label="avg")
    line4 = ax2.plot(gen, param_min, "r--", label="min")
    ax2.set_ylabel(paramName, color="r")
    for tl in ax2.get_yticklabels():
        tl.set_color("r")

    lns = line1 + line2 + line3 + line4 + line5
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="best")
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close(fig)

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

# Plot just the efficiency, fake rate and duplicate for the best individal (HOF) over time
def plotScores(efficiencies, fakeRates, duplicateRates):
    pathlib.Path(plotDirectory).mkdir(parents=True, exist_ok=True)
    plotName = plotDirectory + "Efficiency_FakeRate_DuplicateRate.png"
    # metrics = ["Efficiency", "FakeRate", "DuplicateRate"]
    fig, ax = plt.subplots()
    l1 = ax.plot(efficiencies, "b-", label="efficiency")
    l2 = ax.plot(fakeRates, "r-", label="fakeRate")
    l3 = ax.plot(duplicateRates, "g-", label="duplicateRate")
    ax.axis([0, NGEN - 1, -1, 100])
    lns = l1 + l2 + l3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="best")
    plt.tight_layout()
    # print("Now we should be saving the figure..")
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


def main():
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
    # Objects that will compile the data for population graphs
    logbook = tools.Logbook()
    popData = {}
    popData["Score"] = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    popData["Efficiency"] = tools.Statistics(key=lambda ind: ind.fitness.values[1])
    popData["FakeRate"] = tools.Statistics(key=lambda ind: ind.fitness.values[2])
    popData["DuplicateRate"] = tools.Statistics(key=lambda ind: ind.fitness.values[3])
    for oneName in NAME_TO_FACTOR:
        popData[oneName] = tools.Statistics(key=lambda ind: ind[NAME_TO_INDEX[oneName]])
    # mstats = tools.MultiStatistics(popData)
    # mstats = tools.MultiStatistics(Score=popData["Score"], Efficiency=popData["Efficiency"], FakeRate=popData["FakeRate"], DuplicateRate=popData["DuplicateRate"], sigmaScattering=popData["sigmaScattering"],
    #                                 maxSeedsPerSpM=popData["maxSeedsPerSpM"], maxPt=popData["maxPt"], impactMax=popData["impactMax"], deltaRMax=popData["deltaRMax"], deltaRMin=popData["deltaRMin"], radLengthPerSeed=popData["radLengthPerSeed"])
    stats_score = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    stats_eff = tools.Statistics(key=lambda ind: ind.fitness.values[1])
    stats_fake = tools.Statistics(key=lambda ind: ind.fitness.values[2])
    stats_dup = tools.Statistics(key=lambda ind: ind.fitness.values[3])
    stats_sigmaScattering = tools.Statistics(
        key=lambda ind: ind[NAME_TO_INDEX["sigmaScattering"]])
    stats_maxSeedsPerSPM = tools.Statistics(
        key=lambda ind: ind[NAME_TO_INDEX["maxSeedsPerSpM"]])
    stats_maxPt = tools.Statistics(key=lambda ind: ind[NAME_TO_INDEX["maxPt"]])
    stats_impactMax = tools.Statistics(
        key=lambda ind: ind[NAME_TO_INDEX["impactMax"]])
    stats_deltaRMin = tools.Statistics(key=lambda ind: ind[NAME_TO_INDEX["deltaRMin"]])
    stats_deltaRMax = tools.Statistics(key=lambda ind: ind[NAME_TO_INDEX["deltaRMax"]])
    stats_radLengthPerSeed = tools.Statistics(key=lambda ind: ind[NAME_TO_INDEX["radLengthPerSeed"]])
    mstats = tools.MultiStatistics(Score=stats_score, Efficiency=stats_eff, FakeRate=stats_fake, DuplicateRate=stats_dup, sigmaScattering=stats_sigmaScattering,
                                maxSeedsPerSpM=stats_maxSeedsPerSPM, maxPt=stats_maxPt, impactMax=stats_impactMax, deltaRMax=stats_deltaRMax, deltaRMin=stats_deltaRMin, radLengthPerSeed=stats_radLengthPerSeed)
    mstats.register("avg", np.mean)
    mstats.register("std", np.std)
    mstats.register("min", np.min)
    mstats.register("max", np.max)

    hof = tools.HallOfFame(1)
    # initialize NPOP copies of initial guess
    pop = toolbox.population_guess()
    indPrint(pop[0])
    firstFit = toolbox.evaluate(pop[0])
    print(firstFit)
    for ind in pop:
        ind.fitness.values = firstFit
    g = 0
    bestEff = 0
    bestDup = 20
    bestFake = 100
    # Stop condition is that the efficiency is greater than 99.4, duplicate % is less than 60, and fakerate is less than 10%
    # or the number of generations reaches NGEN (100).
    while g < NGEN and ((bestEff < 99.4) or bestDup > 60 or bestFake > 10):
        g = g + 1
        print("-- Generation %i --" % g)
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
        # Apply mutation on the offspring
        maxMutants = 16
        mutantsCount = 1
        numParamsDiffCount = [] # for debugging
        for i, mutant in enumerate(offspring):
            prevInd = []
            for j in range(len(mutant)):
                prevInd.append(mutant[j])
            if mutantsCount == maxMutants:
                break # I only have 16 cores so evaluating more individuals will slow down
            if random.random() < MUTPB:
                mutantsCount += 1
                toolbox.mutate(mutant)
                numMutatedParams = 0
                for j in range(len(mutant)): # for debugging
                    if prevInd[j] != mutant[j]:
                        numMutatedParams += 1
                del mutant.fitness.values
                numParamsDiffCount.append(numMutatedParams)
        print(f"List of # of mutated params in each individual: {numParamsDiffCount}")
        # Evaluate the individuals with an invalid fitness (the mutated individuals)
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        print(f"Evaluating {len(invalid_ind)} individuals...")
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        # Print and store data
        toPrint = 1
        printCounter = 0
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
            if (printCounter < 1 and fit[1] == -1):
                print("This ind broke the seeding algo")
                indPrint(ind) # print first toPrint individual(s) that cause seedingalgo to break
                printCounter += 1
        
        pop[:] = offspring
        hof.update(pop)
        # Gather all the fitnesses in one list and print the stats
        scoreList = [ind.fitness.values[0] for ind in pop]
        effs = [ind.fitness.values[1] for ind in pop]
        fakeRates = [ind.fitness.values[2] for ind in pop]
        dupRates = [ind.fitness.values[3] for ind in pop]
        length = len(pop)
        mean = sum(effs) / length
        print("Efficiency:")
        printStats(effs, length, "Percent")
        print("Duplicate Rate:")
        printStats(dupRates, length, "Percent")
        print("Fake Rate:")
        printStats(fakeRates, length, "Percent")
        print("The best one so far:", end=" ")
        # record data for analyzing best individual
        goodOne = hof[0]
        bestIndividuals.append(goodOne)
        for oneName in hofData:
            paramVal = goodOne[NAME_TO_INDEX[oneName]] * NAME_TO_FACTOR[oneName]
            if oneName == "maxSeedsPerSpM":
                hofData[oneName].append(int(paramVal))
            else:
                hofData[oneName].append(paramVal)
        scores.append(goodOne.fitness.values[0])
        bestEff = goodOne.fitness.values[1]
        efficiencies.append(bestEff)
        bestFake = goodOne.fitness.values[2]
        fakeRateList.append(bestFake)
        bestDup = goodOne.fitness.values[3]
        dupRateList.append(bestDup)
        indPrint(goodOne)
        print("Best score (Score, efficiency, fakeRate, dupRate):", end=" ")
        print(goodOne.fitness.values)
        # record data for analyzing the population
        logbook.record(gen=g, **mstats.compile(pop))
    # record efficiency over different pT ranges
    for cutRange in pTCutData:
        cutRangeInputs = []
        for goodOne in bestIndividuals:
            names, params = createNamesAndParams(goodOne)
            names.append("fltPrtPtMin")
            params.append(cutRange[0])
            names.append("fltPrtPtMax")
            params.append(cutRange[1])
            args = paramsToInput(params, names)
            cutRangeInputs.append(args)
        print(f"Evaluating pT range {cutRange}")
        avgScores = toolbox.map(executeAlg, cutRangeInputs)
        for avgScore in avgScores:
            pTCutData[cutRange].append(avgScore["efficiency"])
    plotPtRange(pTCutData)
    # Make plots for the population
    for oneName in NAME_TO_FACTOR:
        #print("Length of data is " + str(len(logbook.chapters[oneName].select("max"))))
        plotLogbook("Score", oneName, logbook)
        plotLogbook("Efficiency", oneName, logbook)
        plotLogbook("DuplicateRate", oneName, logbook)
        plotLogbook("FakeRate", oneName, logbook)

    # Make plots for the best individual
    for oneName in hofData:
        plotHOF("Score", scores, oneName, hofData[oneName])
        plotHOF("Efficiency", efficiencies, oneName, hofData[oneName])
        plotHOF("DuplicateRate", dupRateList, oneName, hofData[oneName])
        plotHOF("FakeRate", fakeRateList, oneName, hofData[oneName])
    # Plot the score, efficiency, dupliceate rate and fake rate all here
    plotScores(efficiencies, fakeRateList, dupRateList)
    print(f"Wrote plots to {plotDirectory}")
    return logbook, hof

logbook, hof = main()

print(f"NAME_TO_INDEX = {NAME_TO_INDEX}")
print(f"NAME_TO_FACTOR = {NAME_TO_FACTOR}")
print()
oldGuess = []
for i, oneName in enumerate(NAME_TO_FACTOR):
  oldGuess.append(myGuess[i] * NAME_TO_FACTOR[oneName])
print(f"My guess was {oldGuess}")
if (len(MINS) != len(MAXS)) or len(NAME_TO_FACTOR) != len(MINS):
    print("Mismatched definition of names_to_factor and/or mins maxs")

