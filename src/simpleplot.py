

import matplotlib
from matplotlib.pylab import plt #load plot library
matplotlib.use('Agg')
from sys import stdout
import sys


filepath = str(sys.argv[1])
print('Reading file ', filepath)

plt.figure(1, figsize = (50,20))
axes = plt.gca()
plt.gca().set_aspect('auto')

vector_l = 2000
bin_size = 0
x_size = vector_l
y_size = 500

axes.set_xlim(0, x_size)

word_hash = "none"
frequency = 0
seq_len = 0

count = 0

xdata_prev = 0
ydata_prev = 0

total = 0
observed_prob = 0
p_sum = 0
highest = 0

second_plot = False

with open(filepath) as openfileobject:
	for line in openfileobject:

		# For the second part of the plot (Freqs)
		if(second_plot == True):
			xdata = line.split(",")
			ydata = [highest, highest]
			goodxdata = [int(float(xdata[0])/bin_size), int(float(xdata[1])/bin_size)]
			if(goodxdata[0] == goodxdata[1]):
				goodxdata[1] += 1
			#print(goodxdata)
			plt.plot(goodxdata, ydata, color="purple", linewidth=5)

		if(str(line) == 'EOP\n'):
			second_plot = True
			
		if(count == 0):
			word_hash = line
			print("Desired hash: ", word_hash)

		if(count == 1):
			
			frequency = int(line)
			print('Frequency: ', frequency)
		
		if(count == 2):
			
			seq_len = int(line)
			print('Length of sequence: ', seq_len)
			bin_size = seq_len / vector_l
			print('Bin size: ', bin_size)

		if(count == 3):
			
			observed_prob = float(line)
			print('Observed probability: ', observed_prob)


		if(second_plot == False and count > 3 and count < (vector_l + 3)):

			rvalue = float(line)/float(bin_size)
			if(rvalue > highest):
				highest = rvalue
			p_sum += rvalue

			xdata = [xdata_prev, count-3]
			ydata = [ydata_prev, rvalue]
			plt.plot(xdata, ydata, color="black")

			total = total + int(line)

			#print(ydata)

			xdata_prev = xdata_prev + 1
			ydata_prev = rvalue

		if(second_plot == False and count == (vector_l + 3)):
			xdata_prev = 0
			ydata_prev = 0
		if(second_plot == False and count >= (vector_l + 3)):

			rvalue = float(line)/float(bin_size)

			xdata = [xdata_prev, xdata_prev + 1]
			ydata = [ydata_prev, rvalue]
			plt.plot(xdata, ydata, color="red")

			#print(ydata)

			xdata_prev = xdata_prev + 1
			ydata_prev = rvalue
			
		
			

		count = count + 1
	
print('Average probability density: ', p_sum/vector_l)
plt.title('Frequency concentration of '+word_hash)
plt.xlabel('Coordinates of the sequence')
plt.ylabel('K-mer density [0-1] for word '+word_hash)

avg_freq = float(frequency)/float(seq_len)
uniform = 1 / ( 4 ** len(word_hash) )
print('Average frequency distribution: ', avg_freq)
print('Expected uniform distribution: ', uniform)
print('Total accounted for: ', total)

plt.plot([0, vector_l], [avg_freq, avg_freq], color="blue", label='Average frequency distribution')
plt.plot([0, vector_l], [2*avg_freq, 2*avg_freq], color="cyan", label='Average frequency distribution')
plt.plot([0, vector_l], [uniform, uniform], color="orange", label='Uniform distribution')
plt.plot([0, vector_l], [observed_prob, observed_prob], color="green", label='Observed probability')
leg = axes.legend()


plt.savefig(filepath+'.png')

labels = [item.get_text() for item in axes.get_xticklabels()]
labels = [int(item) * seq_len/vector_l for item in labels]
labels = [ "{:.2E}".format((item))   for item in labels]
axes.set_xticklabels(labels)
plt.tight_layout()

plt.savefig(filepath+'.png')
print('\nSaved figure at ', filepath+'.png')



#plt.plot(x,y)

