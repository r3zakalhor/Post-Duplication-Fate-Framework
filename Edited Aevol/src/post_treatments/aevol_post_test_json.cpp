//
// Created by duazel on 27/05/2020.
//

#include "IOJson.h"
#include "ExpManager.h"
using namespace aevol;

int main(int argc, char **argv) {
  IOJson test("param.in", "chromosome.txt");

  test.load(new ExpManager(), true, NULL, 0, NULL, 0);
  test.write("input.json");


//  IOJson test2("extracted.json");
//  test2.load(new ExpManager(), true, NULL, 0, NULL, 0);
  // test2.write("output.json");
  test.write("extracted.json");


  IOJson test2("extracted.json");
  test2.load(new ExpManager(), true, NULL, 0, NULL, 0);
  test2.write("output.json");
}