#include "unique.h"

#include <vector>
#include <algorithm>
#include <ranges>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

using Catch::Matchers::Equals;

TEST_CASE("Simple") {
    CHECK_THAT(Unique({1, 2, 2, 3, 3, 3}), Equals(std::vector{1, 2, 3}));
    CHECK_THAT(Unique({1, 1, 1, 1}), Equals(std::vector{1}));
    CHECK_THAT(Unique({1, 2, 3, 4}), Equals(std::vector{1, 2, 3, 4}));
    CHECK_THAT(Unique({1, 1, 2, 2, 2, 3}), Equals(std::vector{1, 2, 3}));
    CHECK_THAT(Unique({1, 1}), Equals(std::vector{1}));
}

TEST_CASE("Empty") {
    CHECK(Unique({}).empty());
    CHECK_THAT(Unique({0}), Equals(std::vector{0}));
}

TEST_CASE("Big") {
    auto input = GENERATE(take(10, chunk(100'000, random(-30'000, 30'000))));
    std::ranges::sort(input);

    auto expected = input;
    auto [last, end] = std::ranges::unique(expected);
    expected.erase(last, end);

    CHECK_THAT(Unique(input), Equals(expected));
}
