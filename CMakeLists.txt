cmake_minimum_required(VERSION 3.14)

project(shad-cpp0)

set(CMAKE_CXX_STANDARD             20)
set(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/tools/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS  ON)

find_package(Catch REQUIRED)

include(tools/cmake/TestSolution.cmake)
include(tools/cmake/BuildFlags.cmake)
include_directories(tools/util)

function(add_if_exists name)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name})
    add_subdirectory(${name})
  endif()
endfunction()

add_if_exists(private)

add_if_exists(multiplication)

add_if_exists(bin-pow)
add_if_exists(range)
add_if_exists(unique)
add_if_exists(split)
add_if_exists(polish-notation)
add_if_exists(word-count)
add_if_exists(permutations)
add_if_exists(reverse-map)
add_if_exists(long-sum)
add_if_exists(diff-pairs)

add_if_exists(quadratic)
add_if_exists(pt-triangle)
add_if_exists(factorization)
add_if_exists(rotate)
add_if_exists(sort-students)
add_if_exists(war)
add_if_exists(filter-even)
add_if_exists(swap-sort)
add_if_exists(entrance)

add_if_exists(ring-buffer)
add_if_exists(matrix)
add_if_exists(static-map)
add_if_exists(string-view)
add_if_exists(dungeon)
add_if_exists(stack)

add_if_exists(find-last)
add_if_exists(local-max)
add_if_exists(partition)
add_if_exists(fold)
add_if_exists(functors)
add_if_exists(itertools)

find_package(PNG)
if (${PNG_FOUND})
  add_if_exists(raytracer-debug)
  add_if_exists(raytracer-geom)
  add_if_exists(raytracer-reader)
  add_if_exists(raytracer)
  add_if_exists(raytracer-b2)
elseif()
  message(STATUS "raytracer disabled. PNG:${PNG_FOUND}")
endif()

add_if_exists(advance)
add_if_exists(universal_find)

add_if_exists(vector)
add_if_exists(lru-cache)
add_if_exists(cow-vector)

add_if_exists(dedup)
add_if_exists(smart-ptr)
add_if_exists(deque)

add_if_exists(list)
add_if_exists(calculator)
add_if_exists(reader)

add_if_exists(scheme)
add_if_exists(scheme-tokenizer)
add_if_exists(scheme-parser)

add_if_exists(safe-transform)
add_if_exists(intrusive-list)
add_if_exists(any)

add_if_exists(tryhard)
add_if_exists(pimpl)

find_package(Poco QUIET COMPONENTS Foundation Net JSON)
if (Poco_FOUND)
  add_if_exists(hello-http)
  add_if_exists(bot)
else()
  message(STATUS "Seems like POCO is not installed on your machine.")
  message(STATUS "Tasks that depend on POCO are disabled.")
endif()

add_if_exists(editor)
add_if_exists(small-test-framework)

add_if_exists(scala-vector)
add_if_exists(optional)

find_package(nanoflann QUIET)
if (${PNG_FOUND} AND ${nanoflann_FOUND})
  add_if_exists(fast-dither)
elseif()
  message(STATUS "fast-dither disabled. PNG:${PNG_FOUND} nanoflann:${nanoflann_FOUND}")
endif()

add_if_exists(decoder)
add_if_exists(bad-hash)
add_if_exists(bad-rand)
