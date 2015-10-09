#include <fstream>

#include "catch.hpp"
#include "../../include/Topology.hpp"

TEST_CASE("Parse small Topology from file")
{
	NSMCuda::Topology t = NSMCuda::Topology();
	std::ifstream topology_file("test_data/topology.txt");

	REQUIRE( topology_file.good() );

	REQUIRE_NOTHROW( topology_file >> t );

	SECTION("Check subvolumes number")
	{
		REQUIRE( t.getN() == 4 );
	}

	SECTION("Check neighbours number")
	{
		REQUIRE( t.getNeighboursLength() == 6 );
	}

	SECTION("Check neighbours array")
	{
		int expected[] = {1, 0, 2, 1, 3, 2};

		int nl = t.getNeighboursLength();
		int * arr = t.getNeighboursArray();

		for(int i = 0; i < nl; i++) {
			if (expected[i] != arr[i]) {
				WARN( "failed at index " + std::to_string(i) );
			}
			CHECK( expected[i] == arr[i] );
		}
	}

	SECTION("Check offset array")
	{
		int expected[] = {0, 1, 3, 5};

		int * arr = t.getOffsetArray();

		for(int i = 0; i < t.getN(); i++) {
			if (expected[i] != arr[i]) {
				WARN( "failed at index " + std::to_string(i) );
			}
			CHECK( expected[i] == arr[i] );
		}
	}

}
