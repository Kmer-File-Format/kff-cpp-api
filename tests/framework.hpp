/**
 * Minimal unit test framework (no dependency).
 *
 *   TEST_CASE("group: description") { CHECK(...); CHECK_EQ(a, b); }
 *
 * The runner (main.cpp) executes every registered case, or only the cases
 * whose name contains one of the command line arguments.
 */
#ifndef KFF_TEST_FRAMEWORK_HPP
#define KFF_TEST_FRAMEWORK_HPP

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace kfftest {

struct Failure {
	std::string message;
};

struct Case {
	const char * name;
	void (*function)();
};

inline std::vector<Case> & registry() {
	static std::vector<Case> cases;
	return cases;
}

inline uint64_t & check_counter() {
	static uint64_t counter = 0;
	return counter;
}

struct Registrar {
	Registrar(const char * name, void (*function)()) {
		registry().push_back({name, function});
	}
};

// Printable representation of the compared values
template <typename T>
std::string repr(const T & value) {
	std::ostringstream os;
	os << value;
	return os.str();
}
inline std::string repr(const uint8_t & value) { return std::to_string((unsigned)value); }
inline std::string repr(const int8_t & value) { return std::to_string((int)value); }
inline std::string repr(const char & value) { return std::string("'") + value + "'"; }
inline std::string repr(const std::string & value) { return "\"" + value + "\""; }
inline std::string repr(const std::vector<uint8_t> & bytes) {
	std::ostringstream os;
	os << "[" << bytes.size() << " bytes:";
	for (size_t i = 0; i < bytes.size() && i < 64; i++)
		os << " " << std::hex << std::setw(2) << std::setfill('0') << (unsigned)bytes[i];
	if (bytes.size() > 64)
		os << " ...";
	os << "]";
	return os.str();
}

[[noreturn]] inline void fail(const char * file, int line, const std::string & what) {
	std::ostringstream os;
	os << file << ":" << line << ": " << what;
	throw Failure{os.str()};
}

} // namespace kfftest

#define KFF_TEST_CAT2(a, b) a##b
#define KFF_TEST_CAT(a, b) KFF_TEST_CAT2(a, b)

#define TEST_CASE(name)                                                                   \
	static void KFF_TEST_CAT(kff_test_case_, __LINE__)();                                \
	static kfftest::Registrar KFF_TEST_CAT(kff_test_registrar_, __LINE__)(              \
		name, &KFF_TEST_CAT(kff_test_case_, __LINE__));                                   \
	static void KFF_TEST_CAT(kff_test_case_, __LINE__)()

#define CHECK(condition)                                                                  \
	do {                                                                                  \
		++kfftest::check_counter();                                                       \
		if (!(condition))                                                                 \
			kfftest::fail(__FILE__, __LINE__, "CHECK(" #condition ") failed");            \
	} while (0)

#define CHECK_EQ(actual, expected)                                                        \
	do {                                                                                  \
		++kfftest::check_counter();                                                       \
		const auto & kff_a_ = (actual);                                                   \
		const auto & kff_e_ = (expected);                                                 \
		if (!(kff_a_ == kff_e_))                                                          \
			kfftest::fail(__FILE__, __LINE__,                                             \
				"CHECK_EQ(" #actual ", " #expected ")\n      actual:   " +                  \
				kfftest::repr(kff_a_) + "\n      expected: " + kfftest::repr(kff_e_));      \
	} while (0)

// Checks that the expression throws an exception of the given type
#define CHECK_THROWS_AS(expression, type)                                                 \
	do {                                                                                  \
		++kfftest::check_counter();                                                       \
		bool kff_thrown_ = false;                                                         \
		try { expression; } catch (type const &) { kff_thrown_ = true; }                  \
		if (!kff_thrown_)                                                                 \
			kfftest::fail(__FILE__, __LINE__, "CHECK_THROWS_AS(" #expression ", " #type   \
				") did not throw");                                                      \
	} while (0)

#endif
