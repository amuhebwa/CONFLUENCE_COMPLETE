#!/bin/bash
echo "**************** starting pre-diagnostics ***********************************"
docker run -v /mnt/d/Confluence/CONFLUENCE_COMPLETE:/usr/local/src/app -p 8080:8080 amuhebwa/confluence_prediagnostics
echo "**************** pre-diagnostics completed ***********************************"


echo "**************** starting bam computations ***********************************"
docker run  -v "$PWD":/usr/local/src/app -p 8081:8081 amuhebwa/bam_algorithm

echo "**************** bam computations completed ***********************************"
echo "**************** starting Integrator computations***********************************"
docker run  -v "$PWD":/usr/local/src/app -p 8082:8082 amuhebwa/confluence_integrator
echo "**************** Integrator computations finished ***********************************"

echo "**************** starting post-diagnostics ***********************************"
docker run  -v "$PWD":/usr/local/src/app -p 8083:8083 amuhebwa/confluence_postdiagnostics
echo "**************** post-diagnostics completed ***********************************"

