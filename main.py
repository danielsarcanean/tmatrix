import numpy as np
import matplotlib.pyplot as plt
import csv

###----Constants----###

plt.rcParams["figure.dpi"] = 300
fig, ax = plt.subplots(figsize=(10,5))
me =  1 # SI: 9.11E-31 kg
hbar =  1 / (2 * np.pi) # SI: 6.626E-34 / (2 * np.pi) J · s
energy = 7 / 27.21 # --> convertion to hartrees from eV 

###-----------------###

def tMatrixPositiveEnergy(xa, xb):
    length = positions[xb] - positions[xa]
    
    if energy < potentials[xa]:
        k = np.sqrt(2 * me * (potentials[xa] - energy) / hbar**2)
        matrix = np.array([
            [np.exp(k*length), 0],
            [0, np.exp(k*length)]
        ])
        x = np.linspace(positions[xa], positions[xa] + 10, 200)
        plt.plot(x, (energy)*np.exp(-k*(x)), color='red')
    
    elif energy > potentials[xa] and energy < potentials[xb]:
        k = np.sqrt(2 * me * (potentials[xb] - energy) / hbar**2)
        matrix = np.array([
            [np.exp(k*length), 0],
            [0, np.exp(-k*length)]
        ])
        x = np.linspace(positions[xb], positions[xb] + 10, 200)
        plt.plot(x, (energy)*np.exp(-k*(x - positions[xb])), color='red')
    
    elif potentials[xa] == potentials[xb]:
        k = np.sqrt(2 * me * (energy - potentials[xa]) / hbar**2)
        matrix = np.array([
            [np.exp(1j*k*length), 0],
            [0, np.exp(-1j*k*length)]
        ])
        x = np.linspace(positions[xa], positions[xb], 200)
        #plt.plot(x, energy +(energy - potentials[xb])*np.exp(-1j*k*(x)), color='green')
        plt.plot(x, (energy) + (energy-potentials[xb])*np.exp(1j*k*(x- positions[xb])), color='red')
    
    else:
        k1 = np.sqrt(2 * me * (energy - potentials[xa]) / hbar**2)
        k2 = np.sqrt(2 * me * (energy - potentials[xb]) / hbar**2)
        matrix = (0.5)*np.array([
            [1 + k2/k1, 1 - k2/k1],
            [1 - k2/k1, 1 + k2/k1]
        ])
    
    return matrix

def graphicSolution(potentials, positions):
    for i in range(len(potentials)-1):
        if potentials[i] == potentials[i+1]:
            ax.hlines(y=potentials[i], xmin=positions[i], xmax=positions[i+1])
        else: 
            ax.vlines(x=positions[i], ymin=potentials[i], ymax=potentials[i+1])
    plt.ylabel("Potential (hartrees)")
    plt.xlabel("Position (bohr radius)")
    plt.xlim(min(positions), max(positions))
    plt.ylim(1.5*min(potentials), 1.5*max(potentials))
    plt.legend(frameon=False, fontsize=9)
    plt.show()


def transferParameters(xa, xb):
    finalMatrix = np.array([[1, 0], [0, 1]])
    for i in range(xb - xa):
        nextMatrix = tMatrixPositiveEnergy(i, i+1)
        finalMatrix = np.matmul(finalMatrix, nextMatrix)
        if energy < potentials[i]:
            break
        elif energy > potentials[i] and energy < potentials[i+1]:
            break
        else:
            continue
    transmitance = (finalMatrix[0,0] * finalMatrix[1,1] - finalMatrix[0,1] * finalMatrix [1,0]) / finalMatrix[1,1]
    reflectance = finalMatrix[1,0] / finalMatrix[1,1]
    print(f"Transmitance = {np.abs(transmitance)**2:.3f} | Reflectance = {np.abs(reflectance)**2:.3f}")
    return finalMatrix

### Variables of the program ### 

potentials = []
positions = []

with open('potentials.csv', 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=';', skipinitialspace=True)
        
    for row in csv_reader:
        potentials.append(float(row['Potentials (eV)'])/27.21) # We ensure the eV values are passed to hartrees 
        positions.append(float(row['Positions (nm)'])*18.897) # in bohr radius



### Execution ### 

transferParameters(0,5)
graphicSolution(potentials, positions)
