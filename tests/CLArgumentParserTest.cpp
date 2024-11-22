/*
#include <gtest/gtest.h>


#include <vector>
#include <fstream>

#include "../src/io/CLArgumentParser.h"
#include "../src/forces/LennardJones.h"
#include "../src/io/file/in/CuboidReader.h"
*/
/*
 * @brief Parse example CLI arguments correctly
 */
/*
TEST(CLArgumentParser, parse) {
    std::ofstream file("testFile.txt");
    if (file.is_open()) {
        file << "Hello, World!" << std::endl;
        file.close();
    }

    int argc = 9;

    //this is because of some warning
    char arg0[] = "MolSim";
    char arg1[] = "-f";
    char arg2[] = "testFile.txt";
    char arg3[] = "-t";
    char arg4[] = "1.0";
    char arg5[] = "-d";
    char arg6[] = "1.1";
    char arg7[] = "-l";
    char arg8[] = "warn";

    char* argv[] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8 };

    Arguments arguments = {
        "",
        10,
        0.014,
        1,
        "info",
        std::make_unique<LennardJones>(),
        std::make_unique<CuboidReader>(),
    };

    //this covers -f testFile.txt
    EXPECT_EQ(CLArgumentParser::parse(argc, argv, arguments), 0);

    EXPECT_EQ(arguments.t_end, 1.0);
    EXPECT_EQ(arguments.delta_t, 1.1);
    EXPECT_STREQ(arguments.logLevel.c_str(), "warn");
    
    //TODO expand to cover force selection
}
*/