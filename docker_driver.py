import os
# mycmd1 = os.popen('docker run -v "$PWD":/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics').read()
cmd1 = 'docker run -v "$PWD":/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics'
cmd2 = 'docker run  -v "$PWD":/usr/local/src/app -p 8081:8081 amuhebwa/bam_algorithm'
cmd3 = 'docker run  -v "$PWD":/usr/local/src/app -p 8082:8082 amuhebwa/confluence_integrator'
cmd4 = 'docker run  -v "$PWD":/usr/local/src/app -p 8083:8083 amuhebwa/confluence_postdiagnostics'
os.popen(cmd1).read()
print("=================================================================================")
os.popen(cmd2).read()
print("=================================================================================")
os.popen(cmd3).read()
print("=================================================================================")
os.popen(cmd4).read()
print("***********TADAAAAAAAAAAAAAA***********************")