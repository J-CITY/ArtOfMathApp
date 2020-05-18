#pragma once
#include <map>
#include <string>

struct Lsystem {
	std::string axiom = "";
	std::map<char, std::string> rules;
	int iterations;
	int angle = 90;
	float distance = 0.0;
	float scaleFactor = 1.0f;

	int widthStep = 0;
	int widthLine = 1;
	float angleStep = 0;

	bool swapPlusMinus = false;
	
	std::string calcLsystem() {
		auto start = axiom;

		if (iterations == 0) {
			return axiom;
		}

		std::string end = "";

		for (auto i = 0; i < iterations; i++) {
			end = "";
			for (auto& ch : start) {
				if (rules.find(ch) != rules.end()) {
					end += rules[ch];
				}
				else {
					end += ch;
				}
			}
			start = end;
		}
		return end;
	}
};
