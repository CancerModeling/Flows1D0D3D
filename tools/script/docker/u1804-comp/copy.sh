id=$(docker create u1804-comp)
# docker cp $id:/ - > ./pd.tar.gz
docker cp $id:/build.log - > ./build.log
docker rm -v $id
