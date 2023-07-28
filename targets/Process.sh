
target=$1

cd $target
echo $target
tar -xf free_energy.tar.gz

rm free_energy.tar.gz
cd free_energy

for i in *
do
	cd $i
	if [ ! -f complex.prmtop ]
	then
		echo $(pwd)
	fi
	tar -czf complex.prmtop.tar.gz complex.prmtop
	tar -czf complex.inpcrd.tar.gz complex.inpcrd
	rm complex.prmtop complex.inpcrd
	cd ..
done

cd ..

cd ..

