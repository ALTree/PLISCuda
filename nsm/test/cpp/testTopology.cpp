#include <fstream>

#include "catch.hpp"
#include "../../include/cpp/Topology.hpp"

TEST_CASE("Parse small Topology from file")
{
	NSMCuda::Topology t = NSMCuda::Topology();
	std::ifstream topology_file("test_data/topology.txt");

	REQUIRE(topology_file.good());

	REQUIRE_NOTHROW(topology_file >> t);

	SECTION("Check subvolumes number"){
	REQUIRE( t.getN() == 4 );
}

	SECTION("Check neighbours number"){
	REQUIRE( t.getNeighboursLength() == 24 );
}

	SECTION("Check neighbours array"){
	int expected[] = {
		1, -1, -1, -1, -1, -1,
		0, 2, -1, -1, -1, -1,
		1, 3, -1, -1, -1, -1,
		2, -1, -1, -1, -1, -1,
	};

	int nl = t.getNeighboursLength();
	int * arr = t.getNeighboursArray();

	putchar('\n');

	for(int i = 0; i < nl; i++) {
		if (expected[i] != arr[i]) {
			WARN( "failed at index " + std::to_string(i) );
		}
		CHECK( expected[i] == arr[i] );
	}
}

}
