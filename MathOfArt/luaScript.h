#pragma once
#include <complex>
#include <utility>

#include "lua/src/lua.hpp"
#include "utils.h"

class LuaScript {
    lua_State* luaState = nullptr;
public:
    LuaScript() {
        
    }

	void init() {
        luaState = luaL_newstate();
        static const luaL_Reg lualibs [] = {
            {"base", luaopen_base},
            {"io", luaopen_io},
            {"math", luaopen_math},
            {NULL, NULL}
        };
        const luaL_Reg* lib = lualibs;
        for (; lib->func != NULL; lib++) {
            luaL_requiref(luaState, lib->name, lib->func, 1);
            lua_settop(luaState, 0);
        }
    }

	~LuaScript() {
        //free(luaState);
    }

    std::tuple<std::complex<double>, std::complex<double>, int> runScript(std::string name,
		std::complex<double> z, std::complex<double> c, int maxi) {
        lua_createtable(luaState, 2, 0);
        lua_pushnumber(luaState, 1);
        lua_pushnumber(luaState, z.real());
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 2);
        lua_pushnumber(luaState, z.imag());
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 3);
        lua_pushnumber(luaState, c.real());
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 4);
        lua_pushnumber(luaState, c.imag());
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 5);
        lua_pushnumber(luaState, maxi);
        lua_settable(luaState, -3);
        lua_setglobal(luaState, "arg");

        int status = luaL_loadfile(luaState, name.c_str());

        int result = 0;
        if (status == LUA_OK) {
            result = lua_pcall(luaState, 0, LUA_MULTRET, 0);
        }
        else {
            E_LOG(" Could not load the script");
        }

        int n = lua_tonumber(luaState, lua_gettop(luaState));
        lua_pop(luaState, 1);
        auto ci = lua_tonumber(luaState, lua_gettop(luaState));
        lua_pop(luaState, 1);
        auto cr = lua_tonumber(luaState, lua_gettop(luaState));
        lua_pop(luaState, 1);
        auto zi = lua_tonumber(luaState, lua_gettop(luaState));
        lua_pop(luaState, 1);
        auto zr = lua_tonumber(luaState, lua_gettop(luaState));
        lua_pop(luaState, 1);
        return std::tuple<std::complex<double>, std::complex<double>, int>(
            std::complex<double>(zr, zi),
            std::complex<double>(cr, ci), n);
    }

    void runAttractorScript(std::string name, int maxi, int w, int h, int scale) {
        lua_createtable(luaState, 2, 0);
        lua_pushnumber(luaState, 1);
        lua_pushnumber(luaState, maxi);
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 2);
        lua_pushnumber(luaState, w);
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 3);
        lua_pushnumber(luaState, h);
        lua_settable(luaState, -3);
        lua_pushnumber(luaState, 4);
        lua_pushnumber(luaState, scale);
        lua_settable(luaState, -3);
        lua_setglobal(luaState, "arg");

        int status = luaL_loadfile(luaState, name.c_str());

        int result = 0;
        if (status == LUA_OK) {
            result = lua_pcall(luaState, 0, LUA_MULTRET, 0);
        }
        else {
            E_LOG(" Could not load the script");
        }
    }
};
