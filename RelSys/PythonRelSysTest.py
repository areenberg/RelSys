import relsys

arrivalRates = [0.8,2.5,0.6,2.8]
serviceTimes = [10,5,10,8]
capacity = [15,20,10,30]
relocationProbabilities = [[0.0,0.4,0.1,0.5],
                           [0.3,0.0,0.5,0.0],
                           [0.0,0.5,0.0,0.5],
                           [0.2,0.3,0.5,0.0]]
preferredQueue = [0,1,2,3]

#import parameters
relsys.input(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)
relsys.setType("approximation")
relsys.equalize(True)

#run the model
relsys.run()

#results
#for queueIdx in range(4):
#    print(relsys.getDensity(queueIdx))

#for queueIdx in range(4):
#    print(relsys.getFreq(queueIdx))

for queueIdx in range(4):
    print(relsys.getShortageProb(queueIdx))

#for queueIdx in range(4):
#    print(relsys.getAvailProb(queueIdx))
    
#for queueIdx in range(4):
#    print(relsys.getExpOccupany(queueIdx))

#for queueIdx in range(4):
#    print(relsys.getExpOccFraction(queueIdx))
