
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

  // this is because of some warning
  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "testFile.txt";
  char arg3[] = "-s";
  char arg4[] = "0.5";
  char arg5[] = "-l";
  char arg6[] = "warn";

  char* argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6};  // , arg7, arg8 };

  auto [name, step, write_checkpoint] = CLArgumentParser::parse(argc, argv);
  EXPECT_EQ(name, "testFile.txt");
  EXPECT_EQ(step, 0.5);
  EXPECT_EQ(SpdWrapper::get()->level(), spdlog::level::warn);
  EXPECT_EQ(write_checkpoint, false);

#if 0  // TODO: i dont know why it says file is empty or directoy here even
       // though its the same test
  // checks default write checkpoint
  if (std::ofstream ffile("TestFile.txt"); ffile.is_open()) {
    ffile << "Hello, World!" << std::endl;
    ffile.close();
  }
  char arg7[] = "-c";
  char arg21[] = "testFile.txt";

  char* argv1[] = {arg0, arg1, arg21, arg3, arg4, arg5, arg6}; //, arg7};
  auto [name1, step1, write_checkpoint1] = CLArgumentParser::parse(argc, argv1);
  //EXPECT_EQ(write_checkpoint1, false);
#endif
}
