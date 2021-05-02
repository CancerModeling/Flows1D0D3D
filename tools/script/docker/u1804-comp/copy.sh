id=$(docker create u1804-comp)
docker cp $id:/usr/lib/ - > ./u1804-comp.tar.gz
# docker cp $id:/opt/petsc/3.14.4/configure.log - > ./configure.log
docker rm -v $id
