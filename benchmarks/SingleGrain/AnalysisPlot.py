from pylab import *
import csv

myfile = open("R_2_graph.dat","r")
csv_reader = csv.reader(myfile, delimiter=',', quoting=csv.QUOTE_NONE)
next(csv_reader, None)
step=[]; analytic=[]; simulation=[]

for row in csv_reader:
 step.append(float(row[0]))
 simulation.append(float(row[1]))
 analytic.append(float(row[2]))

fig = plt.figure()
plt.plot(step, analytic,'-', label="Analytic", linewidth = 2)
plt.plot(step, simulation,'-', label="OpenPhase", linewidth = 2)
plt.legend()
plt.ylabel('$R^2$')
plt.xlabel('timesteps')

plt.show()
fig.savefig('R2.png')

myfile.close()
