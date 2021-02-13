# ResolveCompilerPaths - this module defines two macros
#
# RESOLVE_LIBRARIES (XXX_LIBRARIES LINK_LINE)
#  This macro is intended to be used by FindXXX.cmake modules.
#  It parses a compiler link line and resolves all libraries
#  (-lfoo) using the library path contexts (-L/path) in scope.
#  The result in XXX_LIBRARIES is the list of fully resolved libs.
#  Example:
#
#    RESOLVE_LIBRARIES (FOO_LIBRARIES "-L/A -la -L/B -lb -lc -ld")
#
#  will be resolved to
#
#    FOO_LIBRARIES:STRING="/A/liba.so;/B/libb.so;/A/libc.so;/usr/lib/libd.so"
#
#  if the filesystem looks like
#
#    /A:       liba.so         libc.so
#    /B:       liba.so libb.so
#    /usr/lib: liba.so libb.so libc.so libd.so
#
#  and /usr/lib is a system directory.
#
#  Note: If RESOLVE_LIBRARIES() resolves a link line differently from
#  the native linker, there is a bug in this macro (please report it).
#
# RESOLVE_INCLUDES (XXX_INCLUDES INCLUDE_LINE)
#  This macro is intended to be used by FindXXX.cmake modules.
#  It parses a compile line and resolves all includes
#  (-I/path/to/include) to a list of directories.  Other flags are ignored.
#  Example:
#
#    RESOLVE_INCLUDES (FOO_INCLUDES "-I/A -DBAR='\"irrelevant -I/string here\"' -I/B")
#
#  will be resolved to
#
#    FOO_INCLUDES:STRING="/A;/B"
#
#  assuming both directories exist.
#  Note: as currently implemented, the -I/string will be picked up mistakenly (cry, cry)

macro (RESOLVE_LIBRARIES_DIR LIBS_DIR LINK_LINE)
  string (REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" _all_tokens "${LINK_LINE}")
  set (_libs_found)
  set (_directory_list)
  foreach (token ${_all_tokens})
    string (REGEX REPLACE "-L" "" token ${token})
    string (REGEX REPLACE "//" "/" token ${token})
    string (REGEX REPLACE "\n" "" token ${token})
    list (APPEND _directory_list ${token})
  endforeach (token)
  set (${LIBS_DIR} ${_directory_list})
endmacro (RESOLVE_LIBRARIES_DIR)

macro (RESOLVE_LIBRARIES LIBS LINK_LINE)
  #message(${LIBS})
  #message(${LINK_LINE})
  string (REGEX MATCHALL "(((-L|-l|-Wl)([^\" ]+|\"[^\"]+\")|/[^\" ]+(a|so|dll|dylib))|(-framework Accelerate))" _all_tokens "${LINK_LINE}")
  #message(${_all_tokens})
  set (_libs_found)
  set (_directory_list)
  foreach (token ${_all_tokens})
    if (token MATCHES "-L([^\" ]+|\"[^\"]+\")")
      # If it's a library path, add it to the list
      string (REGEX REPLACE "^-L" "" token ${token})
      string (REGEX REPLACE "//" "/" token ${token})
      string (REGEX REPLACE "\n" "" token ${token})
      list (APPEND _directory_list ${token})
    endif (token MATCHES "-L([^\" ]+|\"[^\"]+\")")
  endforeach (token)
  foreach (token ${_all_tokens})
    if (token MATCHES "^(-l([^\" ]+|\"[^\"]+\")|/[^\" ]+(a|so|dll|dylib))")
      string (REGEX REPLACE "\n" "" token ${token})
      # It's a library, resolve the path by looking in the list and then (by default) in system directories
      string (REGEX REPLACE "^-l" "" token ${token})
      set (_root)
      if (token MATCHES "^/")	# We have an absolute path, add root to the search path
	set (_root "/")
      endif (token MATCHES "^/")
      set (_lib "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
      find_library (_lib ${token} HINTS ${_directory_list} ${_root})
      if (_lib)
	      string (REPLACE "//" "/" _lib ${_lib})
        list (APPEND _libs_found ${_lib})
      else (_lib)
        #message (STATUS "Unable to find library ${token}")
      endif (_lib)
    endif (token MATCHES "^(-l([^\" ]+|\"[^\"]+\")|/[^\" ]+(a|so|dll|dylib))")
  endforeach (token)
  foreach (token ${_all_tokens})
    if (token MATCHES "-framework Accelerate")
      list(APPEND _libs_found ${token})
    endif (token MATCHES "-framework Accelerate")
  endforeach (token)
  if (_libs_found)
    list (REVERSE _libs_found)
    list (REMOVE_DUPLICATES _libs_found)
    list (REVERSE _libs_found)
  endif (_libs_found)
  set (${LIBS} "${_libs_found}")
endmacro (RESOLVE_LIBRARIES)

macro (RESOLVE_INCLUDES INCS COMPILE_LINE)
  string (REGEX MATCHALL "-I([^\" ]+|\"[^\"]+\")" _all_tokens "${COMPILE_LINE}")
  set (_incs_found)
  foreach (token ${_all_tokens})
    string (REGEX REPLACE "^-I" "" token ${token})
    string (REGEX REPLACE "//" "/" token ${token})
    string (REGEX REPLACE "\n" "" token ${token})
    if (EXISTS ${token})
      list (APPEND _incs_found ${token})
    else (EXISTS ${token})
      #message (STATUS "Include directory ${token} does not exist")
    endif (EXISTS ${token})
  endforeach (token)
  if (_libs_found)
    list (REMOVE_DUPLICATES _incs_found)
  endif (_libs_found)
  set (${INCS} "${_incs_found}")
endmacro (RESOLVE_INCLUDES)

macro (RESOLVE_DEFINITIONS DEF COMPILE_LINE)
  string (REGEX MATCHALL "-D[0-9A-Z_]+(=[a-zA-z0-9]+)?" _all_tokens "${COMPILE_LINE}")
  # message(STATUS ${_all_tokens})
  set (_incs_found)
  foreach (token ${_all_tokens})
    string (REGEX REPLACE "^-D" "-DMOAB_" token ${token})
    list (APPEND _incs_found ${token})
  endforeach (token)
  if (_libs_found)
    list (REMOVE_DUPLICATES _incs_found)
  endif (_libs_found)
  set (${DEF} "${_incs_found}")
endmacro (RESOLVE_DEFINITIONS)
