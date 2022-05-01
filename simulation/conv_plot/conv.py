import pandas as pd
import matplotlib.pyplot as plt

f = open("./conv_cs50.txt",encoding='utf-8')
plt.style.use('ggplot')
plt.subplots(figsize=[12,4])
ax1 = plt.subplot(1,3,1)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=50")
plt.xlabel("iterations(m)")
plt.ylabel("MSE of $\hat{C}$")

f = open("./conv_cs200.txt",encoding='utf-8')
ax2 = plt.subplot(1,3,2)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=200")
plt.xlabel("iterations(m)")

f = open("./conv_cs500.txt",encoding='utf-8')
ax3 = plt.subplot(1,3,3)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=500")
plt.xlabel("iterations(m)")


plt.tight_layout()
plt.savefig('./ConvCS.png')


f = open("./conv_ar50.txt",encoding='utf-8')
plt.style.use('ggplot')
plt.subplots(figsize=[12,4])
plt.subplot(1,3,1)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=50")
plt.xlabel("iterations(m)")
plt.ylabel("MSE of $\hat{C}$")

f = open("./conv_ar200.txt",encoding='utf-8')
plt.subplot(1,3,2)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=200")
plt.xlabel("iterations(m)")

f = open("./conv_ar500.txt",encoding='utf-8')
plt.subplot(1,3,3)
while True:
    line = f.readline()
    if line:
        line = [eval(i) for i in line.strip().split()]
        plt.plot(line[1:])
    else:
        break
f.close()
plt.title("p=500")
plt.xlabel("iterations(m)")

plt.tight_layout()
plt.savefig('./ConvAR.png')