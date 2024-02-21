# version         v0.1.7
# date            04.02.2024
# author          Alexander-Maurice Illig
# affilation      Institute of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import random
from ._utils import X_to_deltaE
from math import exp
from multiprocessing import Pool

class Explore:
    """
    Description
    -----------
    This class is for exploring the fitness landscape using simulated annealing.
    The exploration conditions can be defined using various parameters (maxSubstitutions, factor, maxIter).
    The search is for maxima; if minima are of interest, the parameter 'sign' should be set to -1.
    
    Attributes
    ----------
    encodeCls : object
        Initialized 'Encode' class.

    model : object
        Initialized and trained 'CombinedPredictor' class.

    yWt : float
        Fitness value of the wild type.
    
    maxSubstitutions : int
        Maximum number of substitutions to introduce (default = 3).

    sign : int
        Either +1 for maxima or -1 for minima (default = +1).

    factor : float
        Accept step if new fitness value is greater than previous*factor (default = 1.0).

    maxIter : int
        Number of iterations to perform (default = 10000).

    """

    def __init__(self, encodeCls:object, model:object, yWt:float, maxSubstitutions=3, sign=+1, factor=1.0, maxIter=1000):
        self.encodeCls = encodeCls
        self.model = model
        self.yWt = yWt
        self.maxSubstitutions = maxSubstitutions
        self.sign = sign
        self.factor = factor
        self.maxIter = maxIter


    def _random_substitution(self) -> str:
        """
        Description
        -----------
        Generate a random substitution to introduce.

        Returns
        -------
        Substitution as string in the format (wildTypeAminoAcid Position variantAminoAcid).
        """
        idx = random.choice(range(self.encodeCls.targetSeq.size))
        wtAA, position = self.encodeCls.targetSeq[idx], self.encodeCls.offsetMap[idx]+self.encodeCls.startingPosition-1
        varAA = random.choice(self.encodeCls.alphabet)
        return '%s%d%s'%(wtAA, position, varAA)

    def _accept_substitution(self, previous, actual, temperature=1e-2) -> bool:
        """
        Description
        -----------
        Metropolis-Hastings algorithm
            1) Introduce random substitution to xt -> x*
            2) Accept if dy = sign*(y* - y*factor) > 0
               Or if r <= min(1, pA)
               with 
               r : uniform random number in [0, 1]
               pA = exp(dy/temperature)
            3) If accept: x* -> xt+1
               Else xt -> xt+1
               
        Parameters
        ---------
        previous : float
            Fitness value of the previous variant.

        actual : float
            Fitness value of the new variant.

        temperature : float
            Temperature to choose (default = 1e-2).
            
        Returns
        -------
            True if substitution should be introduced, False else.
        """
        diff = self.sign*(actual - previous*self.factor)
        
        if diff > 0:
            return True
    
        else:
            if random.uniform(0, 1) <= min(1, exp(diff/temperature)):
                return True

            else:
                return False
    
    @staticmethod
    def _T(step:int, startingTemperature=1e-2, decay=5e-3) -> float:
        """
        Description
        -----------
        Monotonic function for simulated annealing.
        
        Parameters
        ----------
        step : int
            Step in the simulated annealing process.
            
        startingTemperature : float
            Temperature to initialize the process of simulated annealing (default = 1e2).
        
        decay : float
            Strength of 'decay' of monotonic function (default = 5e-3).
        
        Returns
        -------
        Temperature to choose for the 'step' according to the monotonic function.
        """
        return startingTemperature*exp(-decay*step)
    

    def _random_walker(self, seed) -> tuple:
        """
        Description
        -----------
        Generate a "random walker" with random seed 'seed' on the fitness landscape.

        Parameters
        ---------
        seed : int
            Seed to initialize the random generator to create a random walker.

        Returns
        -------
            Tuple including the name and fitness of the improved variant.
            If no improved variant was found, wild type will be returned.
        """
        random.seed(seed)
        
        yPrevious = self.yWt
        
        substitutions = []
        iteration = 0

        while iteration < self.maxIter:
            temperature = self._T(iteration)
            
            iteration += 1

            if len(substitutions) >= self.maxSubstitutions:
                break

            substitution = self._random_substitution()
            variantTemp = substitutions.copy()

            if substitution in substitutions: # exclude same substitution
                continue

            elif substitution[0] == substitution[-1]: # exclude WT
                continue

            elif any([substitution[1:-1] == x[1:-1] for x in variantTemp]): # exclude same position
                variantTemp.remove([x for x in variantTemp if substitution[1:-1] == x[1:-1]][0])

            variantTemp.append(substitution)    

            xVariant = [self.encodeCls._encode_variant(','.join(variantTemp))]
            deltaEVariant = X_to_deltaE(xVariant, self.encodeCls.xWt)
            yActual = self.model.predict(xVariant, deltaEVariant)

            if self._accept_substitution(yPrevious, yActual, temperature):
                substitutions.append(substitution)
                yPrevious = yActual

        if substitutions:
            xVariant = [self.encodeCls._encode_variant(','.join(substitutions))]
            deltaEVariant = X_to_deltaE(xVariant, self.encodeCls.xWt)
            return (','.join(sorted(substitutions, key=lambda x:int(x[1:-1]))), self.model.predict(xVariant, deltaEVariant)[0])
        else:
            return ('WT', self.yWt)

    def scrape_landscape(self, nWalkers, nCores=1) -> set:
        """
        Description
        -----------
        Function to scrape the fitness landscape for improved variants using 'nWalkers'.

        Parameters
        ---------
        nWalkers : int
            Number of random walkers to create to scrape the fitness landscape.

        nCores : int
            Number of cores used for parallelization (default = 1).

        Returns
        -------
            Set including tuples of improved variants and their (predicted) fitness values.
        """
        seeds = range(nWalkers)
        pool = Pool(nCores)

        results = set()
        variants, fitnesses = zip(*pool.map(self._random_walker, seeds))
        for variant, fitness in zip(variants, fitnesses):
            if fitness > self.sign*self.yWt:
                results.add((variant, round(fitness, ndigits=2)))

        return sorted(results, key=lambda x:x[1], reverse=True)
