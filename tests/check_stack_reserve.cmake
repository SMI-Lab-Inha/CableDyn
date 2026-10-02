# SPDX-License-Identifier: Apache-2.0
# Verify the stack reserve written into Windows executables by the GNU link: the 64 MiB
# reserve of the CableDyn executables, and the 1 MiB reserve of the small-stack builds
# that the stack_* tests run (a later --stack option must override the directory-wide one).
#
# CHECKS: '|'-separated <binary>=<expected SizeOfStackReserve in bytes> entries.

foreach(_required OBJDUMP CHECKS)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "check_stack_reserve: missing ${_required}")
  endif()
endforeach()

string(REPLACE "|" ";" _checks "${CHECKS}")
set(_failures "")
foreach(_check IN LISTS _checks)
  string(FIND "${_check}" "=" _eq REVERSE)
  if(_eq LESS 1)
    message(FATAL_ERROR "check_stack_reserve: malformed entry '${_check}'")
  endif()
  string(SUBSTRING "${_check}" 0 ${_eq} _binary)
  math(EXPR _value_at "${_eq} + 1")
  string(SUBSTRING "${_check}" ${_value_at} -1 _expected)
  if(NOT EXISTS "${_binary}")
    message(FATAL_ERROR "check_stack_reserve: missing ${_binary}")
  endif()
  execute_process(
    COMMAND "${OBJDUMP}" -p "${_binary}"
    RESULT_VARIABLE _rc
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "check_stack_reserve: ${OBJDUMP} -p ${_binary} failed (${_rc})\n${_err}")
  endif()
  if(NOT _out MATCHES "SizeOfStackReserve[ \t]+([0-9a-fA-F]+)")
    message(FATAL_ERROR "check_stack_reserve: no SizeOfStackReserve read from ${_binary}")
  endif()
  math(EXPR _reserve "0x${CMAKE_MATCH_1}")
  get_filename_component(_name "${_binary}" NAME)
  if(NOT _reserve EQUAL _expected)
    string(APPEND _failures "  ${_name}: stack reserve ${_reserve} bytes, expected ${_expected}\n")
  else()
    message(STATUS "${_name}: stack reserve ${_reserve} bytes")
  endif()
endforeach()

if(_failures)
  message(FATAL_ERROR "check_stack_reserve: unexpected stack reserve\n${_failures}")
endif()
message(STATUS "PASS: every checked executable has the expected stack reserve")
