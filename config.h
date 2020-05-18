#pragma once
#pragma warning(disable : 4996) 
#include <string>
#include "3rd/cpptoml.h"
#include "utils.h"
#include <chrono>
#include <ctime>

#include "lsystem.h"

struct Config {
	enum class Type {
		MANDELBROT,
		JULIA,
		NEWTON,
		LSYSTEM,
		LAPUNOV
	};

	//main
	Type type = Type::MANDELBROT;
	int width = 1920;
	int height = 1080;
	std::string path = "";
	std::string outFilename = "";
	bool isVideo = false;
	bool isSlides = false;
	unsigned int threadSize = 8;
	int copyFrame = 1;
	int fps = 1;
	//M, J
	unsigned int maxIteration = 212;
	int dim = 2;
	CenterPoint centerStart, centerEnd;
	CenterPoint range = CenterPoint(1.0, 1.0);
	CenterPoint rangeTo = CenterPoint(1.0, 1.0);
	CenterPoint rangeStep = CenterPoint(0.0, 0.0);
	CenterPoint centerStep;
	rgb bgColor;
	bool interactiveMode = false;
	std::string script = "";

	double rotateAngleStart = 0.0;
	double rotateAngleEnd = 0.0;
	double rotateAngleStep = 0.0;

	bool useCinJulia = true;
	double r = 0.0;
	double a = 0.0;
	double aTo = 0.0;
	double aStep = 0.0;
	std::complex<double> c;

	bool isCloud = false;
	int cloudStep = 1;
	//Newton
	std::vector<double> coefs;
	int polynomDegree = 0;

	//lapunov
	double xn = 0.5;
	std::string lapunovPattern = "A";
	
	//Lsystem
	int iterations = 1;
	int iterationsTo = 1;
	int distance = 1;
	double angle = 0.0;
	double scaleFactor = 0.0;
	std::string axiom = "";
	std::vector<std::string> rules;
	bool isStepByStep = false;
	int lineWidth = 1;
	int lineWidthStep = 0;
	int angleStep = 0;
	
	//colors
	std::string colorImage = "";
	std::vector<std::vector<rgb>> colorsMat;
	int minColorRange = 0;
	int maxColorRange = 700;
	int countStepColorX = 10;
	int countStepColorY = 10;
	int stepSizeX = 1;
	int stepSizeY = 1;
	
	Config(std::string _path = "") : path(_path) {
		if (path != "") {
			readConfig(path);
		}
	}
	
	void readConfig(std::string path) {
		try {
			auto config = cpptoml::parse_file(path);
			//[main]
			width = readParam<int>(config, "main.width", width);
			height = readParam<int>(config, "main.height", height);
			auto t = readParam<std::string>(config, "main.type", "MANDELBROT");
			type = strToType[t];
			isVideo = readParam<bool>(config, "main.is_video", isVideo);
			isSlides = readParam<bool>(config, "main.is_slides", isSlides);
			interactiveMode = readParam<bool>(config, "main.interactive_mode", interactiveMode);
			outFilename = readParam<std::string>(config, "main.out_path", outFilename);
			if (outFilename == "") {
				auto t = std::chrono::system_clock::now();
				std::time_t _time = std::chrono::system_clock::to_time_t(t);
				outFilename = std::ctime(&_time);
			}
			threadSize = readParam<int>(config, "main.thread_size", threadSize);
			copyFrame = readParam<int>(config, "main.frames_copy", copyFrame);
			fps = readParam<int>(config, "main.fps", fps);

			auto mainTable = config->get_table("main");
			auto vals = mainTable->get_array_of<int64_t>("bg_color");
			if (vals && (*vals).size() == 3) {
				bgColor = rgb((*vals)[0], (*vals)[1], (*vals)[2]);
			}
			else {
				I_LOG("Can't find key: main.bg_color");
			}

			if (type == Config::Type::MANDELBROT) {
				runMorJParams(config, "mandelbrot");
				isCloud = readParam<bool>(config, "mandelbrot.is_cloud", isCloud);
				cloudStep = readParam<int>(config, "mandelbrot.cloud_step", cloudStep);
				rotateAngleStart = readParam<double>(config, "mandelbrot.rotate_angle", rotateAngleStart);
				rotateAngleEnd = readParam<double>(config, "mandelbrot.rotate_to", rotateAngleEnd);
				rotateAngleStep = readParam<double>(config, "mandelbrot.rotate_step", rotateAngleStep);
				script = readParam<std::string>(config, "mandelbrot.script", script);
			}
			else if (type == Config::Type::JULIA) {
				runMorJParams(config, "julia");
				r = readParam<double>(config, "julia.r", r);
				a = readParam<double>(config, "julia.a", a);
				aTo = readParam<double>(config, "julia.a_to", aTo);
				aStep = readParam<double>(config, "julia.a_step", aStep);

				isCloud = readParam<bool>(config, "julia.is_cloud", isCloud);
				cloudStep = readParam<int>(config, "julia.cloud_step", cloudStep);

				rotateAngleStart = readParam<double>(config, "julia.rotate_angle", rotateAngleStart);
				rotateAngleEnd = readParam<double>(config, "julia.rotate_to", rotateAngleEnd);
				rotateAngleStep = readParam<double>(config, "julia.rotate_step", rotateAngleStep);

				script = readParam<std::string>(config, "julia.script", script);
				
				auto table = config->get_table("julia");
				auto valsC = table->get_array_of<double>("c");
				if (valsC && valsC->size() == 2) {
					c.real((*valsC)[0]);
					c.imag((*valsC)[1]);
				}
				else {
					useCinJulia = false;
					I_LOG("Can't find key: julia.c (will be use r)");
				}
			}
			else if (type == Config::Type::LAPUNOV) {
				runMorJParams(config, "lapunov");
				xn = readParam<double>(config, "lapunov.xn", xn);
				lapunovPattern = readParam<std::string>(config, "lapunov.pattern", lapunovPattern);
			}
			else if (type == Config::Type::NEWTON) {
				runMorJParams(config, "newton");
				auto table = config->get_table("newton");
				auto valsC = table->get_array_of<double>("coefs");
				if (valsC) {
					for (auto& c : *valsC) {
						coefs.push_back(c);
					}
					polynomDegree = coefs.size() - 1;
				}
				else {
					useCinJulia = false;
					I_LOG("Can't find key: newton.coefs");
				}
				rotateAngleStart = readParam<double>(config, "newton.rotate_angle", rotateAngleStart);
				rotateAngleEnd = readParam<double>(config, "newton.rotate_to", rotateAngleEnd);
				rotateAngleStep = readParam<double>(config, "newton.rotate_step", rotateAngleStep);

			}
			else if (type == Config::Type::LSYSTEM) {
				auto lsconf = readParam<std::string>(config, "lsystem.rule", "");
				if (lsconf != "") {
					parseLsystem(lsconf);
				}
			}

			//[color]
			auto clrconf = readParam<std::string>(config, "color.sheme", "");
			if (clrconf != "") {
				parseColor(clrconf);
			}
		}
		catch (const cpptoml::parse_exception& e) {
			E_LOG("Error parse: " + path);
		}
	}

	Lsystem getLsystem() {
		Lsystem ls;
		ls.iterations = iterations;
		ls.widthLine = lineWidth;
		ls.widthStep = lineWidthStep;
		ls.angle = angle;
		ls.angleStep = angleStep;
		ls.axiom = axiom;
		ls.distance = distance;
		ls.scaleFactor = scaleFactor;
		for (auto& r : rules) {
			auto rule = split(r, ':');
			ls.rules[rule[0][0]] = rule[1];
		}
		
		return ls;
	}

private:
	template<typename T>
	T readParam(std::shared_ptr<cpptoml::table>& config, std::string paramName, T defVal) {
		if (auto param = config->get_qualified_as<T>(paramName)) {
			return param.value_or(defVal);
		}
		I_LOG("Can't find key: " + paramName);
		return defVal;
	}

	void runMorJParams(std::shared_ptr<cpptoml::table>& config, std::string t) {
		dim = readParam<int>(config, t+".dim", dim);
		maxIteration = readParam<int>(config, t + ".max_iterations", maxIteration);

		auto table = config->get_table(t);

		auto valsRange = table->get_array_of<double>("range");
		if (valsRange && valsRange->size() == 2) {
			range.x = (*valsRange)[0];
			range.y = (*valsRange)[1];
		}
		else {
			I_LOG("Can't find key: " + t + ".range");
		}
		auto valsRangeTo = table->get_array_of<double>("range_to");
		if (valsRangeTo && valsRangeTo->size() == 2) {
			rangeTo.x = (*valsRangeTo)[0];
			rangeTo.y = (*valsRangeTo)[1];
		}
		else {
			I_LOG("Can't find key: " + t + ".range_to");
		}
		auto valsRangeStep = table->get_array_of<double>("range_step");
		if (valsRangeStep && valsRangeStep->size() == 2) {
			rangeStep.x = (*valsRangeStep)[0];
			rangeStep.y = (*valsRangeStep)[1];
		}
		else {
			I_LOG("Can't find key: " + t + ".range_step");
		}
		
		auto vals = table->get_array_of<double>("center");
		if (vals && vals->size() == 2) {
			centerStart.x = (*vals)[0];
			centerStart.y = (*vals)[1];
		}
		else {
			I_LOG("Can't find key: " + t + ".center");
		}
		auto valsTo = table->get_array_of<double>("center_to");
		if (valsTo && valsTo->size() == 2) {
			centerEnd.x = (*valsTo)[0];
			centerEnd.y = (*valsTo)[1];
		}
		else {
			I_LOG("Can't find key: " + t + ".center_to");
		}
		auto valsStep = table->get_array_of<double>("center_step");
		if (valsStep && valsStep->size() == 2) {
			centerStep.x = (*valsStep)[0];
			centerStep.y = (*valsStep)[1];
		}
		else {
			I_LOG("Can't find key: "+ t + ".center_step");
		}
	}
	
	void parseLsystem(std::string lspath) {
		try {
			auto config = cpptoml::parse_file(lspath);

			iterations = readParam<int>(config, "iterations", iterations);
			iterationsTo = readParam<int>(config, "iterations_to", iterations);
			distance = readParam<int>(config, "distance", distance);
			angle = readParam<double>(config, "angle", angle);
			scaleFactor = readParam<double>(config, "scale_factor", scaleFactor);
			axiom = readParam<std::string>(config, "axiom", axiom);
			isStepByStep = readParam<bool>(config, "is_step_by_step", isStepByStep);

			lineWidth = readParam<int>(config, "line_width", lineWidth);
			lineWidthStep = readParam<int>(config, "line_width_step", lineWidthStep);
			angleStep = readParam<int>(config, "angle_step", angleStep);
			
			auto valsRules = config->get_array_of<std::string>("rules");
			if (valsRules) {
				for (auto v : *valsRules) {
					rules.push_back(v);
				}
			}
			else {
				I_LOG("Can't find key: " + lspath + " : rules");
			}
		}
		catch (const cpptoml::parse_exception& e) {
			E_LOG("Error parse: " + path);
		}
	}
	
	void parseColor(std::string cpath) {
		try {
			auto config = cpptoml::parse_file(cpath);

			colorImage = readParam<std::string>(config, "image", colorImage);
			if (colorImage != "") {
				return;
			}
			
			//minColorRange = readParam<int>(config, "min_range", minColorRange);
			//maxColorRange = readParam<int>(config, "max_range", maxColorRange);
			countStepColorX = readParam<int>(config, "count_step_x", countStepColorX);
			countStepColorY = readParam<int>(config, "count_step_y", countStepColorY);
			stepSizeX = readParam<int>(config, "step_size_x", stepSizeX);
			stepSizeY = readParam<int>(config, "step_size_y", stepSizeY);
			
			auto valsColors = config->get_array_of<cpptoml::array>("colors");
			if (valsColors) {
				for (auto v : *valsColors) {
					colorsMat.push_back(std::vector<rgb>());
					auto cs = v->get_array_of<cpptoml::array>();
					if (!cs) {
						E_LOG("Wrong color array");
						continue;
					}
					for (auto c : *cs) {
						auto rgbarr = c->get_array_of<int64_t>();
						if (!rgbarr) {
							E_LOG("Wrong color array");
							continue;
						}
						if (rgbarr->size() != 3) {
							colorsMat[colorsMat.size() - 1].push_back(rgb());
							E_LOG("Wrong color array size (!= 3)");
							continue;
						}
						colorsMat[colorsMat.size() - 1].push_back(rgb((*rgbarr)[0], (*rgbarr)[1], (*rgbarr)[2]));
					}
				}
			}
			else {
				I_LOG("Can't find key: " + cpath + " : colors");
			}
		}
		catch (const cpptoml::parse_exception& e) {
			E_LOG("Error parse: " + path);
		}
	}

	std::map<std::string, Type> strToType {
		{"MANDELBROT", Type::MANDELBROT},
		{"JULIA", Type::JULIA},
		{"NEWTON", Type::NEWTON},
		{"LSYSTEM", Type::LSYSTEM},
		{"LAPUNOV", Type::LAPUNOV}
	};
};
