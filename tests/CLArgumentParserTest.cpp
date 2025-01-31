
#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include "../src/io/CLArgumentParser.h"

/*
 * @brief Parse example CLI arguments correctly
 */
TEST(CLArgumentParser, parse) {
  if (std::ofstream file("testFile.txt"); file.is_open()) {
    file << "Hello, World!" << std::endl;
    file.close();
  }

  int argc = 7;

  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "testFile.txt";
  char arg3[] = "-s";
  char arg4[] = "0.5";
  char arg5[] = "-l";
  char arg6[] = "warn";

  char* argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6};

  auto [name, step, checkpoint_path] = CLArgumentParser::parse(argc, argv);
  EXPECT_EQ(name, "testFile.txt");
  EXPECT_EQ(step, 0.5);
  EXPECT_EQ(SpdWrapper::get()->level(), spdlog::level::warn);
  EXPECT_FALSE(checkpoint_path);
}
