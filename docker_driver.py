import subprocess
import time
import code
# for unix based environment. It seems split function messes us the path interpreted as $PWD
#cmd1 = 'docker run -v "$PWD":/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics'
#cmd2 = 'docker run  -v "$PWD":/usr/local/src/app -p 8081:8081 amuhebwa/bam_algorithm'
#cmd3 = 'docker run  -v "$PWD":/usr/local/src/app -p 8082:8082 amuhebwa/confluence_integrator'
#cmd4 = 'docker run  -v "$PWD":/usr/local/src/app -p 8083:8083 amuhebwa/confluence_postdiagnostics'
#cmd1 = 'docker run -v /Users/amuhebwa/Documents/Python/CONFLUENCE_COMPLETE:/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics'
#cmd2 = 'docker run -v /Users/amuhebwa/Documents/Python/CONFLUENCE_COMPLETE:/usr/local/src/app -p 8081:8081 amuhebwa/bam_algorithm'
#cmd3 = 'docker run -v /Users/amuhebwa/Documents/Python/CONFLUENCE_COMPLETE:/usr/local/src/app -p 8082:8082 amuhebwa/confluence_integrator'
#cmd4 = 'docker run -v /Users/amuhebwa/Documents/Python/CONFLUENCE_COMPLETE:/usr/local/src/app -p 8083:8083 amuhebwa/confluence_postdiagnostics'

# For windows
cmd1 = 'docker run -v C:\CONFLUENCE_COMPLETE:/usr/local/src/app amuhebwa/confluence_prediagnostics'
cmd2 = 'docker run -v C:\CONFLUENCE_COMPLETE:/usr/local/src/app amuhebwa/bam_algorithm'
cmd3 = 'docker run -v C:\CONFLUENCE_COMPLETE:/usr/local/src/app amuhebwa/confluence_integrator'
cmd4 = 'docker run -v C:\CONFLUENCE_COMPLETE:/usr/local/src/app amuhebwa/confluence_postdiagnostics'


print("************* STARTING PRE-DIAGNOSTICS ***************")
process1 = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
output1, error1 = process1.communicate()
print(process1)
del process1, output1, error1
print("************* PRE-DIAGNOSTICS COMPLETED ***************")
time.sleep(5)

print("************* STARTING BAM COMPUTATIONS ***************")

process2 = subprocess.Popen(cmd2.split(), stdout=subprocess.PIPE)
output2, error2 = process2.communicate()

print(process2)
del process2, output2, error2
print("************* BAM COMPUTATIONS COMPLETED ***************")
time.sleep(5)
print("************* STARTING INTEGRATOR CALCULATIONS **************")
process3 = subprocess.Popen(cmd3.split(), stdout=subprocess.PIPE)
output3, error3 = process3.communicate()

print(process3)
del process3, output3, error3
print("************* INTEGRATOR CALCULATIONS COMPLETED ***************")
time.sleep(5)
print("************* STARTING POST-DIAGNOSTICS ***************")
process4 = subprocess.Popen(cmd4.split(), stdout=subprocess.PIPE)
output4, error4 = process4.communicate()

print(process4)
del process4, output4, error4
print("************* POST-DIAGNOSTICS COMPLETED ***************")
