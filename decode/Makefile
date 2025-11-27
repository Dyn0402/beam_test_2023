all:
	 g++ decode.cpp -o decode $(shell root-config --cflags) $(shell root-config --libs)
	 g++ convert_vec_tree_to_array.cpp -o convert_vec_tree_to_array $(shell root-config --cflags) $(shell root-config --libs)
	 g++ convert_array_tree_to_vec.cpp -o convert_array_tree_to_vec $(shell root-config --cflags) $(shell root-config --libs)
	 g++ subtract_pedestals.cpp -o subtract_pedestals $(shell root-config --cflags) $(shell root-config --libs)

clean:
	 rm -f decode convert_vec_tree_to_array convert_array_tree_to_vec subtract_pedestals
