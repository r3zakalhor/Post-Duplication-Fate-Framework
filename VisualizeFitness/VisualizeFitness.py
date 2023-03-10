# Importing required libraries
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 16
# Creating a series of data of in range of 1-50.
fig, ax = plt.subplots(figsize=(12, 5))
x = np.linspace(0, 1, 200)


# Creating a Function.
def normal_dist(x, h, mean, sd):
    prob_density = h * np.exp(-0.5 * ((x - mean) / sd) ** 2)
    return prob_density

# Calculate mean and Standard deviation.
mean = np.mean(x)
sd = np.std(x)
pdf = 0
pdf2 = 0
pdf3 = 0
pdf4 = 0

l_list = []

if len(sys.argv) < 3:
    print("VisualizeFitnees.py firstparam.in color1 secondparam.in color2")
# Apply function to the data.
with open(sys.argv[1], 'r') as file:
    # reading each line
    for line in file:
        # reading each word
        for word in line.split():
            # displaying the words
            if word == "ENV_ADD_GAUSSIAN":
                # print(word)
                for word1 in line.split():
                    l_list.append(word1)
                # print(l_list)
                pdf += normal_dist(x, float(l_list[1]), float(l_list[2]), float(l_list[3]))
            l_list.clear()

if len(sys.argv) > 3:
    with open(sys.argv[3], 'r') as file:
        # reading each line
        for line in file:
            # reading each word
            for word in line.split():
                # displaying the words
                if word == "ENV_ADD_GAUSSIAN":
                    # print(word)
                    for word1 in line.split():
                        l_list.append(word1)
                    # print(l_list)
                    pdf2 += normal_dist(x, float(l_list[1]), float(l_list[2]), float(l_list[3]))
                l_list.clear()

if len(sys.argv) > 5:
    with open(sys.argv[5], 'r') as file:
        # reading each line
        for line in file:
            # reading each word
            for word in line.split():
                # displaying the words
                if word == "ENV_ADD_GAUSSIAN":
                    # print(word)
                    for word1 in line.split():
                        l_list.append(word1)
                    # print(l_list)
                    pdf3 += normal_dist(x, float(l_list[1]), float(l_list[2]), float(l_list[3]))
                l_list.clear()

    with open(sys.argv[7], 'r') as file:
        # reading each line
        for line in file:
            # reading each word
            for word in line.split():
                # displaying the words
                if word == "ENV_ADD_GAUSSIAN":
                    # print(word)
                    for word1 in line.split():
                        l_list.append(word1)
                    # print(l_list)
                    pdf4 += normal_dist(x, float(l_list[1]), float(l_list[2]), float(l_list[3]))
                l_list.clear()



# pdf = normal_dist(x, 1.2, 0.52, 0.12) + normal_dist(x, -1.4, 0.5, 0.07) + normal_dist(x, 0.3, 0.8, 0.03)
# pdf2 = normal_dist(x, 1.1, 0.62, 0.22) + normal_dist(x, -1.5, 0.6, 0.17) + normal_dist(x, 0.2, 0.9, 0.13)
# Plotting the Results


if len(sys.argv) > 5:
    ax.plot(x, pdf, color=sys.argv[2], label='Env. (a)')
    ax.plot(x, pdf2, color=sys.argv[4], label='Env. (b)')
    ax.plot(x, pdf3, color=sys.argv[6], label='Env. (c)')
    ax.plot(x, pdf4, color=sys.argv[8], label='Env. (d)')

elif len(sys.argv) > 3:
    ax.plot(x, pdf, color=sys.argv[2], label='Env. (a)')
    ax.plot(x, pdf2, color=sys.argv[4], label='Env. (b)')

else:
    ax.plot(x, pdf, color=sys.argv[2], label='Env. (a)')
    
ax.set_ylim(0, 1)
plt.xlabel('Celluar')
plt.ylabel('Efficacy')
plt.legend()


plt.savefig("Fitness curves.pdf", dpi=1000)
print("done.")
plt.show()
