from relsys import relsys

class model:
    
    def __init__(self,
                 arrivalRates,
                 serviceTimes,
                 capacity,
                 relocationProbabilities,
                 preferredQueue):
        
        self.mdl = relsys.model()        
        self.mdl.input(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)
        
    def run(self):
        self.mdl.run()
        
    def getDensity(self,queueIdx=0,type="preferred"):
        return self.mdl.getDensity(queueIdx,type=type)
        
    def getShortageProb(self,queueIdx=0):
        return self.mdl.getShortageProb(queueIdx)
    
    def setType(self,mdltype):
        self.mdl.setType(mdltype)
    
    def queuesEval(self,qEvalIdx):
        self.mdl.queuesEval(qEvalIdx)
    
    def queuesEval(self,qEvalIdx):
        self.mdl.queuesEval(qEvalIdx)
        
    def equalize(self,equalize):
        self.mdl.equalize(equalize)
    
    def setVerbose(self,set):
        self.mdl.setVerbose(set)
    
    def setSeed(self,sd):
        self.mdl.setSeed(sd)
    
    def setAccSamType(self,stype):
        self.mdl.setAccSamType(stype)
    
    def setSimTolerance(self,tol):
        self.mdl.setSimTolerance(tol)
    
    def setBurnIn(self,bin):
        self.mdl.setBurnIn(bin)
    
    def setSimTime(self,mnTime):
        self.mdl.setSimTime(mnTime)
        
    def setSamples(self,mnSamples):
        self.mdl.setSamples(mnSamples)
        
    def setHyperPhases(self,openStates,blockedStates):
        self.mdl.setHyperPhases(openStates,blockedStates)
        
    def getFreq(self,queueIndex=0,type="preferred"):
        return self.mdl.getFreq(queueIndex,type=type)
        
    def getAvailProb(self,queueIndex=0,type="preferred"):
        return self.mdl.getAvailProb(queueIndex,type=type)
    
    def getExpOccupany(self,queueIndex=0):
        return self.mdl.getExpOccupany(queueIndex)
        
    def getExpOccFraction(self,queueIndex=0):
        return self.mdl.getExpOccFraction(queueIndex)
        
    def getArrivalRates(self):
        return self.mdl.getArrivalRates()
    
    def getServiceTimes(self):
        return self.mdl.getServiceTimes()
    
    def getCapacity(self):
        return self.mdl.getCapacity()
        
    def getReloc(self):
        return self.mdl.getReloc()
        
    def getPreferredQueue(self):
        return self.mdl.getPreferredQueue()