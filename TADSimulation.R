#Script Created by Kanishk Asthana on 30 Oct 2014
#Creating New Empty Interaction Matrix with Baseline Interaction Frequency of 30

interaction_matrix<- matrix(rpois(10000,lambda=40), nrow=100,ncol=100)

#Setting Self Interaction of a specific element with itself as Zero

    for(i in 1:100)
    {
        interaction_matrix[i,i]=0
    }

#Printing a Sample
#print(interaction_matrix[1:10,1:10])

#Generating Lengths of TADs using a Poission Distribution such that Total Length of Elements in TADs is always 100

sumoftads=0
tad_lengths=numeric(length=4)

    while(sumoftads>100 || sumoftads==0)
    {    
        #Generating TAD lengths using a Poisson Model
        nums=rpois(3,lambda=25)
        sumoftads=sum(nums)
 #       print(nums)
        tad_lengths[1:3]=nums
    }

tad_lengths[4]=100-sumoftads

print("Actual TAD boundaries:")
print(cumsum(tad_lengths[1:3]))

#TAD Lengths

#print(tad_lengths)

start_pos=1
end_pos=0
    
    for(i in 1:length(tad_lengths))
    {
        end_pos=end_pos + tad_lengths[i];
        for(j in start_pos:end_pos)
        {
            interaction_matrix[j,start_pos:end_pos]=rpois(length(start_pos:end_pos),lambda=90)    
        }
        start_pos=end_pos+1
    }

#Setting Self Interaction of a specific element with itself as Zero

for(i in 1:100)
{
    interaction_matrix[i,i]=0
}

interaction_heatmap <- heatmap(interaction_matrix, Rowv=NA, Colv=NA, col = rev(heat.colors(100)),scale="none",revC=TRUE, main="Simulated HeatMap of Interactions")

print('Computing bias scores first Run')
#Bias scores are computed as follows:
# The sum of all elements upstream of a certain position is stored in a variable and the same is done for the downstream
end=nrow(interaction_matrix)
bias1=numeric(length=end)

    for(i in 1:end)
        {
          upstream_interaction_count=sum(interaction_matrix[i,1:i])
          downstream_interaction_count=sum(interaction_matrix[i,i:end])
          bias1[i]=upstream_interaction_count-downstream_interaction_count
        }

plot(1:end,bias1,type='S',xlab='Position',ylab='Bias Score',col='blue',main='Computing Bias score by summing all the values of all upstream and downstream\n interactions scores and taking the difference of those summed scores')
lines(1:end,numeric(length=end),col='red')
plot(2:end-1,diff(bias1),type='l',ylab='Rate of Change of Bias',xlab='Position',main="Simulated TAD boundaries can be found out by mapping\n rate of change of bias score")

print('Positions of Start/End of TADs assuming they are continuous using Bias score:')
print((2:end-1)[diff(bias1)< -200])

#Although taking the rate of change can be used to find the TAD boundaries, the bias scores do not accurately account for bias. In the ideal case the bias should be
# zero in the middle of the TAD, which is something we do not see. Moreover, the rate of change approach might not not be very useful given greater
# changes in the rate of the TAD boundaries. 

# It may be more appropriate to change the Bias scoring so that interactions scores along only a window rather than the whole chromosome are summed.

bias2=numeric(length=end)
#defining a window size
window=10

    for(i in 1:end)
        {
         if(i<=window)
             {
              upstream_interaction_count=sum(interaction_matrix[i,1:i])
              downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
              bias2[i]=upstream_interaction_count-downstream_interaction_count
             }
         
         else if(i+window>=end)
             {
              upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
              downstream_interaction_count=sum(interaction_matrix[i,i:end])
              bias2[i]=upstream_interaction_count-downstream_interaction_count
             }
         else{
             upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
             downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
             bias2[i]=upstream_interaction_count-downstream_interaction_count
            }
        }
#Plotting new Bias scores
plot(1:end,bias2,type='s',xlab="Position",ylab='Bias Score',main='Summing in a Window leads to better Bias score:\n Here the Bias near the Center of the TAD is close to Zero\nHowever there are some end-effects')
lines(1:end,numeric(length=end),col='red')
#Testing output
#print(bias2)
plot(2:end-1,diff(bias2),type='l',ylab='Rate of Change of Windowed Bias',xlab='Position',main="Simulated TAD boundaries can be found out by mapping\n rate of change of new bias score as well")

print('Positions of Start/End of TADs assuming they are continuous using new Bias score:')
print((2:end-1)[diff(bias2)< -200])

# Calculating the Bias Score using by calculating the ratios of the upstream and downstream interaction scores

bias3=numeric(length=end)

for(i in 1:end)
{
    upstream_interaction_count=sum(interaction_matrix[i,1:i])
    downstream_interaction_count=sum(interaction_matrix[i,i:end])
    bias3[i]=upstream_interaction_count/downstream_interaction_count
}

plot(1:end,bias3,type='S',xlab='Position',ylab='Bias Score',col='blue',main='Computing Bias score by summing all the values of all upstream and downstream\n interactions scores and taking the ratio of those summed scores')
lines(1:end,numeric(length=end),col='red')

#The Bias score does not explain bias very clearly, moreover taking ratio does not seem to be a good way findout where the boundaries are as this is 
#not very clear from the graph, and may pose problems when there is more background noise or other confounding factors
plot(2:end-1,diff(bias3),type='l',ylab='Rate of Change of Bias',xlab='Position',main="Simulated TAD boundaries can be found out by mapping\n rate of change of ratio bias score but less reliably")
print('Position found out using Ratio bias Scoring :')
print((2:end-1)[diff(bias3)< -0.05])


#Like the difference case, the total ratio measure shows end effects and not good bias measurements, where there should be a zero bias or a bias of 1 in the
#middle of the chromosomes

bias4=numeric(length=end)
#defining a window size
window=10

for(i in 1:end)
{
    if(i<=window)
    {
        upstream_interaction_count=sum(interaction_matrix[i,1:i])
        downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
        bias4[i]=upstream_interaction_count/downstream_interaction_count
    }
    
    else if(i+window>=end)
    {
        upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
        downstream_interaction_count=sum(interaction_matrix[i,i:end])
        bias4[i]=upstream_interaction_count/downstream_interaction_count
    }
    else{
        upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
        downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
        bias4[i]=upstream_interaction_count/downstream_interaction_count
    }
}
#Plotting new Bias scores
plot(1:end,bias4,type='s',xlab="Position",ylab='Bias Score',main='Summing in a Window leads to better Bias score:\n Here the Bias near the Center of the TAD is close to One\nHowever there are some end-effects')
lines(1:end,numeric(length=end),col='red')

#Testing output
#print(bias2)
plot(2:end-1,diff(bias4),type='l',ylab='Rate of Change of Windowed Bias',xlab='Position',main="Simulated TAD boundaries can be found out by mapping\n rate of change of new ratio bias score as well")

print("The boundaries are more clearly defined with windowing")
print('Position found out using windowed Ratio bias Scoring :')
print((2:end-1)[diff(bias4)< -1])

#Testing results for Log of ratio:

plot(1:end,log(bias4),type='s',xlab='Position',ylab='Log of Windowed Ration Bias Score',main='Taking log of the ratio of the windowed Bias score makes the center of the TAD \nhave a score near zero and provides better resolution')
lines(1:end,numeric(length=end),col='red')
plot(2:end-1,diff(log(bias4)),type='l',ylab='Rate of change of log of windowed ratio bias score',xlab='Position',
     main='Simulated TAD boundaries can also be found out with this approach quite reliably')

print('Position found out using log windowed Ratio bias Scoring :')
print((2:end-1)[diff(log(bias4))< -1])


