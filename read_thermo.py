import matplotlib.pyplot as plt
import os

if not os.path.exists('thermo'):
    os.makedirs('thermo')

system_name = 'PVP5 Ag(111) NPT'

timestep = []
poteng = []
toteng = []
temp = []
vol = []
press = []

with open('thermo.lammps') as file:
    count = 0
    for line in file:
        if count!=0:
            i = map(float, line.split())
            timestep.append(i[0])
            poteng.append(i[1])
            toteng.append(i[3])
            temp.append(i[4])
            vol.append(i[5])
            press.append(i[6])
        count = count + 1

plt.plot(timestep, poteng, label='PotEng')
plt.plot(timestep, toteng, label='TotEng')
plt.legend()
plt.ylabel('Energy (eV)')
plt.xlabel('Timestep')
plt.title(system_name)
plt.savefig('thermo/energy.png', bbox_inches='tight')
plt.show()
plt.close()

plt.plot(timestep, temp)
plt.ylabel('Temperature (K)')
plt.xlabel('Timestep')
plt.title(system_name)
plt.savefig('thermo/temperature.png', bbox_inches='tight')
plt.show()
plt.close()

plt.plot(timestep, vol)
plt.ylabel(r'Volume ($\AA^3$)')
plt.xlabel('Timestep')
plt.title(system_name)
plt.savefig('thermo/volume.png', bbox_inches='tight')
plt.show()
plt.close()

plt.plot(timestep, press)
plt.ylabel('Pressure (atm)')
plt.xlabel('Timestep')
plt.title(system_name)
plt.savefig('thermo/pressure.png', bbox_inches='tight')
plt.show()
plt.close()
            