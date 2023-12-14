#include <iostream>
#include <lua.hpp>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void TestLua(CCTK_ARGUMENTS)
{
    // create new Lua state
    lua_State *lua_state;
    lua_state = luaL_newstate();

    // load Lua libraries
    static const luaL_Reg lualibs[] =
    {
        { "base", luaopen_base },
        { NULL, NULL}
    };

    const luaL_Reg *lib = lualibs;
    for(; lib->func != NULL; lib++)
    {
        lib->func(lua_state);
        lua_settop(lua_state, 0);
    }

    // run the Lua script
    int status = luaL_loadstring(lua_state, "print('Hello World from Lua!'); return 'abc',42;");

    if (status != LUA_OK)
    {
      CCTK_ERROR("Loading Lua function failed");
    }
    lua_pcall(lua_state, 0, LUA_MULTRET, 0);

    if (!lua_gettop(lua_state))
    {
      CCTK_ERROR("Lua function return didn't work");
    }
    
    // Arguments are reverse-ordered, so start with last
    if (lua_type(lua_state, lua_gettop(lua_state)) != LUA_TNUMBER)
    {
      CCTK_ERROR("Lua function return had wrong type");
    }
    if (lua_tonumber(lua_state, lua_gettop(lua_state)) != 42)
    {
      CCTK_ERROR("Lua function return didn't return correct value.");
    }
    lua_pop(lua_state, 1);

    // Next argument (previous argument in lua-return)
    if (lua_type(lua_state, lua_gettop(lua_state)) != LUA_TSTRING)
    {
      CCTK_ERROR("Lua function return had wrong type");
    }
    {
      std::string tmp (lua_tostring(lua_state, lua_gettop(lua_state)));
      if (tmp.compare("abc"))
      {
        CCTK_ERROR("Lua function return didn't return correct value.");
      }
    }
    CCTK_INFO("Lua return worked.");

    // close the Lua state
    lua_close(lua_state);
}

