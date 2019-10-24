echo 1 > lock_spack

. share/spack/setup-env.sh

mkdir -p mofem_mirror &&
curl -s -L https://bitbucket.org/likask/mofem-cephas/downloads/mirror_v0.9.0.tar.gz \
| tar xzC $PWD/mofem_mirror  --strip 1

spack mirror remove mofem_mirror  || true
spack mirror add mofem_mirror $PWD/mofem_mirror
spack bootstrap
spack install --only dependencies mofem-cephas+slepc

rm -f lock_spack

#if [ ! $GIT_COMMIT == $GIT_PREVIOUS_COMMIT ]; then

	rm -rf /var/lib/jenkins/workspace/SpackBuildDevelopBranch/build
  rm -rf /var/lib/jenkins/workspace/SpackBuildUserModulesDevelopBranch/build

#fi