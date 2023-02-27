#load the relsys module
import relsys 

#----------------------------------
#   MODEL DATA
#----------------------------------

#arrival rates of each customer type
arrivalRates = [0.8,2.5,0.6,2.8]

#mean service time of each customer type
serviceTimes = [10,5,10,8]

#capacity of each queue
capacity = [15,20,10,30]

#fraction of rejected customers that are moved to an alternative queue node
#this is a number of customers x number of queues matrix
relocationProbabilities = [[0.0,0.4,0.1,0.5],
                           [0.3,0.0,0.5,0.0],
                           [0.0,0.5,0.0,0.5],
                           [0.2,0.3,0.5,0.0]]

#queue indices preferred by each customer type
preferredQueue = [0,1,2,3]

#----------------------------------
#   MODEL EVALUATION
#----------------------------------

#import the parameters
relsys.input(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)

#run the model
relsys.run()

#----------------------------------
#   GET THE RESULTS
#----------------------------------

#check the resulting occupancy distribution of each queue 
for queueIdx in range(4):
    print(relsys.getDensity(queueIdx))

#check the resulting shortage probabilities of each queue 
for queueIdx in range(4):
    print(relsys.getShortageProb(queueIdx))