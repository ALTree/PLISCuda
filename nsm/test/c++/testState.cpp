#include <fstream>

#include "catch.hpp"
#include "../../include/State.hpp"

TEST_CASE("Parse small System State from file")
{
	NSMCuda::State s = NSMCuda::State();
	std::ifstream state_file("test_data/state.txt");

	REQUIRE( state_file.good() );

	REQUIRE_NOTHROW( state_file >> s );

	SECTION("Check subvolumes number")
	{
		REQUIRE( s.getN() == 5 );
	}

	SECTION("Check species number")
	{
		REQUIRE( s.getS() == 4 );
	}

	SECTION("Check state array")
	{
		int expected[] = {
				0, 3, 5, 3, 5,
				0, 5, 3, 5, 0,
				5, 0, 0, 0, 0,
				3, 0, 0, 0, 3,
		};

		int size = s.getN() * s.getS();
		int * arr = s.getState();

		for(int i = 0; i < size; i++) {
			if (expected[i] != arr[i]) {
				WARN( "failed at index " + std::to_string(i) );
			}
			CHECK( expected[i] == arr[i] );
		}
	}


}
