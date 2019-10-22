import os
import code
cmd1 = 'docker run -v "$PWD":/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics'
cmd2 = 'docker run  -v "$PWD":/usr/local/src/app -p 8081:8081 amuhebwa/bam_algorithm'
cmd3 = 'docker run  -v "$PWD":/usr/local/src/app -p 8082:8082 amuhebwa/confluence_integrator'
cmd4 = 'docker run  -v "$PWD":/usr/local/src/app -p 8083:8083 amuhebwa/confluence_postdiagnostics'

#list_of_commands = [cmd1, cmd2, cmd3, cmd4]
#code.interact(local=locals())
print("************* STARTING PRE-DIAGNOSTICS ***************")
os.popen(cmd1).read()
print("************* PRE-DIAGNOSTICS COMPLETED ***************")

print("************* STARTING BAM COMPUTATIONS ***************")
os.popen(cmd2).read()
print("************* BAM COMPUTATIONS COMPLETED ***************")

print("************* STARTING INTEGRATOR CALCULATIONS **************")
os.popen(cmd3).read()

print("************* INTEGRATOR CALCULATIONS COMPLETED ***************")

print("************* STARTING POST-DIAGNOSTICS ***************")
os.popen(cmd4).read()
print("************* POST-DIAGNOSTICS COMPLETED ***************")

