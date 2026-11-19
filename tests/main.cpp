#include <cstring>
#include <exception>
#include <iostream>
#include <string>

#include "framework.hpp"

/**
 * Runs all the test cases, or the ones whose name contains one of the
 * command line arguments. Returns the number of failed cases.
 */
int main(int argc, char * argv[]) {
	int failed = 0;
	int run = 0;

	for (const kfftest::Case & test : kfftest::registry()) {
		bool selected = argc <= 1;
		for (int i = 1; i < argc && !selected; i++)
			selected = std::strstr(test.name, argv[i]) != nullptr;
		if (!selected)
			continue;

		run += 1;
		try {
			test.function();
		} catch (const kfftest::Failure & failure) {
			failed += 1;
			std::cerr << "FAILED: " << test.name << "\n  " << failure.message << std::endl;
		} catch (const std::exception & e) {
			failed += 1;
			std::cerr << "FAILED: " << test.name << "\n  unexpected exception: " << e.what() << std::endl;
		} catch (const char * message) {
			failed += 1;
			std::cerr << "FAILED: " << test.name << "\n  unexpected exception: " << message << std::endl;
		} catch (...) {
			failed += 1;
			std::cerr << "FAILED: " << test.name << "\n  unexpected exception" << std::endl;
		}
	}

	std::cout << run - failed << "/" << run << " test cases passed ("
	          << kfftest::check_counter() << " checks)" << std::endl;
	return failed;
}
