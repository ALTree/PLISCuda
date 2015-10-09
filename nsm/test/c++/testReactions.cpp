#include <fstream>

#include "catch.hpp"
#include "../../include/Reactions.hpp"

TEST_CASE("Parse small Reactions system from file")
{
	NSMCuda::Reactions rs = NSMCuda::Reactions();
	std::ifstream reactions_file("test_data/reactions.txt");

	REQUIRE( reactions_file.good() );

	REQUIRE_NOTHROW( reactions_file >> rs );

	SECTION("Check reactions number")
	{
		REQUIRE( rs.getR() == 5 );
	}

	SECTION("Check species number")
	{
		REQUIRE( rs.getS() == 4 );
	}

	SECTION("Check reactants array")
	{
		int expected[] = {
				0, 0, 0, 1, 0,
				0, 0, 1, 0, 1,
				0, 1, 0, 0, 1,
				1, 0, 0, 0, 0,
		};

		int size = rs.getR() * rs.getS();
		int * arr = rs.getReactants();

		for(int i = 0; i < size; i++) {
			if (expected[i] != arr[i]) {
				WARN( "failed at index " + std::to_string(i) );
			}
			CHECK( expected[i] == arr[i] );
		}
	}

	SECTION("Check products array")
	{
		int expected[] = {
				0, 0, 1, 0, 1,
				0, 1, 0, 0, 1,
				1, 0, 0, 0, 1,
				0, 0, 0, 1, 1,
		};

		int size = rs.getR() * rs.getS();
		int * arr = rs.getProducts();

		for(int i = 0; i < size; i++) {
			if (expected[i] != arr[i]) {
				WARN( "failed at index " + std::to_string(i) );
			}
			CHECK( expected[i] == arr[i] );
		}
	}

}
