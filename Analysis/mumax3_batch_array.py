#********************************************************************
#
#  Description:
#
#   This python script generates mumax3 input scripts for a batch
#   simulation on the VCMA skyrmion gate. It uses the mumax3 template
#   script "homework2_natch_template.txt" and substitutes the
#   placeholders with actual values.
#
#   This file is written for the mumax3 workshop organised for the
#   spintronic seminar series in the summer of 2020. For more
#   information go to https://www.mumax.ugent.be/mumax3-workshop
#
#  Authors:
#
#   Dr. Jonathan Leliaert (jonathan.leliaert@ugent.be)
#   Dr. Jeroen Mulkers    (jeroen.mulkers@ugent.be)
#
#********************************************************************
"""
    Description 2:

        Modified to use for sweeping loop steps for
        faster turnaround in hpc3

        NOTE:   template mx3 must have $VARNAME placeholders
                matching that of substitute dict below
                -> don't forget to drop the loop brackets!
"""
from string import Template
import numpy as np
import os


# specify the values for which we want to generate the scripts
# specify filename, this script must be in the same place as template

### COMMENTS FOR RUN ###
comment = "Doing high field run, keeping usual parameters, check whether saturation is correct."

### UNITS OF ergs/cm^2
#start_ku1 = 10800
#end_ku1 = 10800
#steps_ku1 = 1

#start_DMI = 2e-5
#end_DMI = 7e-5
#steps_DMI = 4

#start_msat = 180e3
#end_msat = 195e3
#steps_msat = 4

#1## UNITS OF m
#start_width = 10e-6
#end_width = 12e-6
#steps_width = 3

start_j1 = -5e9
end_j1 = -11e9
steps_j1 = 4

start_j2 = -3e12
end_j2 = -9e12
steps_j2 = 4

#start_Aex = 6e-12
#end_Aex = 9e-12
#steps_Aex = 2

### UNITS OF mA
#start_idc = 0.5
#end_idc = 3.0
#steps_idc = 20

#ku1_values = np.linspace(start_ku1, end_ku1, steps_ku1)
#DMI_values = np.linspace(start_DMI, end_DMI, steps_DMI)
#msat_values = np.linspace(start_msat, end_msat, steps_msat)
j1_values = np.linspace(start_j1, end_j1, steps_j1)
#j1_values = np.linspace(start_j1, end_j1, steps_j1)
j2_values = np.linspace(start_j2, end_j2, steps_j2)
#Aex_values = np.linspace(start_Aex, end_Aex, steps_Aex)
#idc_values = np.linspace(start_idc, end_idc, steps_idc)
#print('DMI ARRAY:\n', DMI_values)
print('RUN COUNT:\t', len(j2_values)*len(j1_values))

template_file = "YIG_SHO_BATCH.mx3"
base_filename = "YIG_SHO"
# read template script
with open(template_file,'r') as f:
    scripttmpl = Template(f.read())

run_num = 1
for j2 in j2_values:
	for j1 in j1_values:
		# write the script for each combination of values
		script = scripttmpl.substitute(dict(J2=j2, J1=j1, COMM=comment))
		scriptfile = base_filename + "j2_%.4f_j1_%.4f__run%g.mx3" % (j2, j1, run_num)
		with open(scriptfile,'w') as f:
		    f.write(script)
		print(scriptfile)
		run_num += 1

print('ARRAY OF SCRIPTS HAVE BEEN MADE')
print('YOU MUST UPDATE THE TASK NUMBER OF RUNS IN THE SBATCH FILE')
