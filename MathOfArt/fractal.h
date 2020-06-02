#pragma once
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <random>
#include <windows.h>
#include "colorManager.h"
#include "config.h"
#include "drawMachine.h"
#include "luaScript.h"

class Fractal {
	cv::VideoWriter* video = nullptr;
public:
	Fractal() {}
	Fractal(std::string cfg) {
		config.readConfig(cfg);
	}
	~Fractal() {
		if (config.isVideo) {
			video->release();
			delete video;
		}
	}

	void save(int id) {
		if (video) {
			for (auto i = 0; i < config.copyFrame; i++) {
				video->write(image);
			}
		}
		else {
			cv::imwrite(config.outFilename + std::to_string(id) + ".png", image);
		}
	}

	int iterations = 0;

	
	void show() {
		centerX = config.centerStart.x;
		centerY = config.centerStart.y;
		range = config.range;
		r = config.r;
		a = config.a;
		initImage();
		//RunJulia();
		//RunMandelbrot();
		//RunLsystem();
		cv::imshow("Display Window", image);

		//imshow("Display Window1", dm.image);
		cv::waitKey(0);
	}

	void run() {
		colorManager = ColorManager(config.seed);
		if (config.isVideo) {
			video = new cv::VideoWriter(config.outFilename+".mp4",
				cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), config.fps,
				cv::Size(config.width, config.height));
		}
		initImage();
		if (config.colorImage != "") {
			colorManager.genColorMatFromImage(config.colorImage);
		}
		else {
			colorManager.genColorGradientMat(config.colorsMat);
		}
		colorManager.resizeColorGrad(config.width, config.height);
		if (config.type == Config::Type::MANDELBROT) {
			auto i = 0;
			centerX = config.centerStart.x;
			centerY = config.centerStart.y;
			range = config.range;
			while (true) {
				if (centerX < config.centerEnd.x)
					centerX = config.centerStart.x + i * config.centerStep.x;
				if (centerY < config.centerEnd.y)
					centerY = config.centerStart.y + i * config.centerStep.y;
				if (range.x < config.rangeTo.x)
					range.x = config.range.x + i * config.rangeStep.x;
				if (range.y < config.rangeTo.y)
					range.y = config.range.y + i * config.rangeStep.y;
				if (angle < config.rotateAngleEnd)
					angle = config.rotateAngleStart + i * config.rotateAngleStep;
				RunMandelbrot();
				save(i);
				i++;

				if ((centerX >= config.centerEnd.x || config.centerStep.x == 0) &&
					(centerY >= config.centerEnd.y || config.centerStep.y == 0) &&
					(range.x >= config.rangeTo.x || config.rangeStep.x == 0) && 
					(range.y >= config.rangeTo.y || config.rangeStep.y == 0) &&
					(angle >= config.rotateAngleEnd || config.rotateAngleStep == 0)) {
					break;
				}
			}
		}
		else if (config.type == Config::Type::JULIA) {
			auto i = 0;
			centerX = config.centerStart.x;
			centerY = config.centerStart.y;
			range = config.range;
			r = config.r;
			a = config.a;
			while (true) {
				if (centerX < config.centerEnd.x)
					centerX = config.centerStart.x + i * config.centerStep.x;
				if (centerY < config.centerEnd.y)
					centerY = config.centerStart.y + i * config.centerStep.y;
				if (range.x < config.rangeTo.x)
					range.x = config.range.x + i * config.rangeStep.x;
				if (range.y < config.rangeTo.y)
					range.y = config.range.y + i * config.rangeStep.y;
				if (a < config.aTo)
					a = config.a + i * config.aStep;
				if (angle < config.rotateAngleEnd)
					angle = config.rotateAngleStart + i * config.rotateAngleStep;
				RunJulia();
				save(i);
				i++;

				if ((centerX >= config.centerEnd.x || config.centerStep.x == 0) && 
					(centerY >= config.centerEnd.y || config.centerStep.y == 0) &&
					(range.x >= config.rangeTo.x || config.rangeStep.x == 0) &&
					(range.y >= config.rangeTo.y || config.rangeStep.y == 0) &&
					(a >= config.aTo || config.aStep == 0) &&
					(angle >= config.rotateAngleEnd || config.rotateAngleStep == 0)) {
					break;
				}
			}
		}
		else if (config.type == Config::Type::NEWTON) {
			auto i = 0;
			centerX = config.centerStart.x;
			centerY = config.centerStart.y;
			range = config.range;
			while (true) {
				if (centerX < config.centerEnd.x)
					centerX = config.centerStart.x + i * config.centerStep.x;
				if (centerY < config.centerEnd.y)
					centerY = config.centerStart.y + i * config.centerStep.y;
				if (range.x < config.rangeTo.x)
					range.x = config.range.x + i * config.rangeStep.x;
				if (range.y < config.rangeTo.y)
					range.y = config.range.y + i * config.rangeStep.y;
				if (angle < config.rotateAngleEnd)
					angle = config.rotateAngleStart + i * config.rotateAngleStep;
				RunFractalNewton();
				save(i);
				i++;

				if ((centerX >= config.centerEnd.x || config.centerStep.x == 0) &&
					(centerY >= config.centerEnd.y || config.centerStep.y == 0) &&
					(range.x >= config.rangeTo.x || config.rangeStep.x == 0) &&
					(range.y >= config.rangeTo.y || config.rangeStep.y == 0) &&
					(angle >= config.rotateAngleEnd || config.rotateAngleStep == 0)) {
					break;
				}
			}
		}
		else if (config.type == Config::Type::LAPUNOV) {
			auto i = 0;
			centerX = config.centerStart.x;
			centerY = config.centerStart.y;
			range = config.range;
			while (true) {
				if (centerX < config.centerEnd.x)
					centerX = config.centerStart.x + i * config.centerStep.x;
				if (centerY < config.centerEnd.y)
					centerY = config.centerStart.y + i * config.centerStep.y;
				if (range.x < config.rangeTo.x)
					range.x = config.range.x + i * config.rangeStep.x;
				if (range.y < config.rangeTo.y)
					range.y = config.range.y + i * config.rangeStep.y;
				RunLapunov();
				save(i);
				i++;

				if ((centerX >= config.centerEnd.x || config.centerStep.x == 0) &&
					(centerY >= config.centerEnd.y || config.centerStep.y == 0) &&
					(range.x >= config.rangeTo.x || config.rangeStep.x == 0) &&
					(range.y >= config.rangeTo.y || config.rangeStep.y == 0)) {
					break;
				}
			}
		}
		else if (config.type == Config::Type::LSYSTEM) {
			auto i = 0;
			iterations = config.iterations;
			while (true) {
				if (iterations < config.iterationsTo)
					iterations++;
				RunLsystem();
				save(i);
				i++;

				if (iterations >= config.iterationsTo) {
					break;
				}
			}
		}
		else if (config.type == Config::Type::BROWNIAN) {
			Brownian();
			save(0);
		}
		else if (config.type == Config::Type::DLA) {
			DLA();
			save(0);
		}
		else if (config.type == Config::Type::ATTRACTOR) {
			Attractor();
			save(0);
		}
		else if (config.type == Config::Type::LIGHTNING) {
			Lightning();
			save(0);
		}
		else if (config.type == Config::Type::LANDSCAPE) {
			if (config.is2d) {
				Landscape2d();
			}
			else {
				Landscape3d();
			}
			save(0);
		}
		else if (config.type == Config::Type::FERN) {
			Fern();
			save(0);
		}
	}
private:
	cv::Mat image;
	Config config;

	void initImage() {
		image = cv::Mat(config.height, config.width, CV_8UC3, 
			cv::Scalar(config.bgColor.blue, config.bgColor.green, config.bgColor.red));
		// cv::Mat::zeros(config.height, config.width, CV_8UC3);
	}

	double angle = 0.0;
	double maxX, minX, minY, centerX = 0, centerY = 0, dx, dy, maxY;
	CenterPoint range = CenterPoint(4.0, 4.0);
	int threadsDone = 0;

	ColorManager colorManager;
	void RunMandelbrot() {
		threadsDone = 0;
		maxX = centerX + range.x / 2;
		minX = centerX - range.x / 2;
		minY = centerY - range.y / 2;
		maxY = centerY + range.y / 2;
		dx = range.x / config.width;
		dy = range.y / config.height;

		int h = (int)(config.height / config.threadSize);
		int dh = config.height - h * config.threadSize;
		for (auto i = 0; i < config.threadSize; ++i) {
			std::thread t(&Fractal::RunMandelbrotThread, this, h * i, (i == config.threadSize - 1 ? h * (i + 1) + dh : h * (i + 1)));
			t.detach();
		}
		while (threadsDone != config.threadSize) {
			Sleep(10);
		}
	}

	double module_sqr(std::complex<double> in)const {
		return in.real() * in.real() + in.imag() * in.imag();
	}

	void RunMandelbrotThread(int sj, int ej) {
		LuaScript luaScript;
		if (config.script != "") {
			luaScript.init();
		}
		double x, y;
		register int n;
		double c = cos(angle), s = sin(angle);
		rgb claudColor;
		double scaleX = (config.width - 1) / range.x;
		double scaleY = (config.height - 1) / range.y;
		for (int j = sj; j < ej; j += config.cloudStep) {
			y = minY + j * dy;
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i += config.cloudStep) {
				x = maxX - i * dx;
				std::complex<double> z((x * c - y * s), (x * s + y * c));
				std::complex<double> res(0, 0);
				n = 0;
				
				claudColor = colorManager.getColorByPos(i, j);

				if (config.script != "" && !config.isCloud) {
					auto zc = luaScript.runScript(config.script, res, z, config.maxIteration);
					res = std::get<0>(zc);
					z = std::get<1>(zc);
					n = std::get<2>(zc);
				}
				else {
					while (module_sqr(res) < 4 && n < config.maxIteration) {
						if (config.script == "") {
							res = _calcMondelbrot(res, config.dim) + z;
						}
						else {
							auto zc = luaScript.runScript(config.script, res, z, config.maxIteration);
							res = std::get<0>(zc);
							z = std::get<1>(zc);
						}

						if (config.isCloud) {
							int _x = res.real() * scaleX + (-(minX + range.x / 2) + config.width / 2);
							int _y = res.imag() * scaleY + (-(minY + range.y / 2) + config.height / 2);

							if (n < 10 || _x < 0 || _y < 0 || _x >= config.width || _y >= config.height) {
								n++;
								continue;
							}

							image.at<cv::Vec3b>(cv::Point(_x, _y)) =
								cv::Vec3b(claudColor.blue, claudColor.green, claudColor.red);
						}

						n++;
					}
				}
				if (!config.isCloud) {
					if (n >= config.maxIteration) {
						pixels[i][2] = 0;
						pixels[i][1] = 0;
						pixels[i][0] = 0;
					}
					else {
						double smoothed = log2(log2(res.real() * res.real() + res.imag() * res.imag()) / 2);
						int colorI = (int)(sqrt(n + 10 - smoothed) * 256) % image.cols;
						//auto color = colorManager.getRgb(colorI);
						auto color = colorManager.getColorByPos(colorI < 0 ? 0 : colorI, 0);
						//auto color = defaultColorSheme(n);
						pixels[i][2] = color.red;
						pixels[i][1] = color.green;
						pixels[i][0] = color.blue;
					}
				}
			}
		}
		threadsDone++;
	}

	rgb defaultColorSheme(int n) {
		double t = (double)n / (double)config.maxIteration;
		int r = (int)(9 * (1 - t) * t * t * t * 255);
		int g = (int)(15 * (1 - t) * (1 - t) * t * t * 255);
		int b = (int)(8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);
		return rgb(r, g, b);
	}

	double _clamp(double s, double e, double val) {
		return val < s ? s : (val > e ? e : val);
	}

	std::complex<double> _calcMondelbrot(std::complex<double>& in, int dim) {
		for (auto i = 1; i < dim; i++) {
			in *= in;
		}
		return in;
	}

	std::complex<double> c;
	double r = 0;
	double a = 0;

	void RunJulia() {
		threadsDone = 0;

		c = config.c;
		if (config.useCinJulia) {
			r = getR(c);
		}
		else {
			c = std::complex<double>(r * cos(a * M_PI / 180.0), 
				r * sin(a * M_PI / 180.0));
		}

		maxX = centerX + range.x / 2;
		minX = centerX - range.x / 2;
		minY = centerY - range.y / 2;
		maxY = centerY + range.y / 2;

		dx = fabs(maxX - minX) / config.width;
		dy = fabs(maxY - minY) / config.height;

		int h = (int)(config.height / config.threadSize);
		int dh = config.height - h * config.threadSize;
		for (auto i = 0; i < config.threadSize; ++i) {
			std::thread t(&Fractal::RunJuliaThread, this, h * i,
				(i == config.threadSize - 1 ? h * (i + 1) + dh : h * (i + 1)));
			t.detach();
		}
		while (threadsDone != config.threadSize) {
			Sleep(10);
		}
	}

	double getR(std::complex<double> c) {
		return (1 + sqrt(1 + 4 * mod(c))) / 2;
	}

	double mod(std::complex<double> c) {
		return sqrt(std::real(c) * std::real(c) + std::imag(c) * std::imag(c));
	}

	void RunJuliaThread(int sj, int ej) {
		LuaScript luaScript;
		if (config.script != "") {
			luaScript.init();
		}
		int n = 0;
		double _c = cos(angle), _s = sin(angle);

		rgb claudColor;
		double scaleX = (config.width - 1) / range.x;
		double scaleY = (config.height - 1) / range.y;
		for (int j = sj; j < ej; j+=config.cloudStep) {
			double y = minY + j * dy;
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i+=config.cloudStep) {
				double x = maxX - i * dx;
				std::complex<double> z((x * _c - y * _s), (x * _s + y * _c));
				int kid = config.maxIteration;

				claudColor = colorManager.getColorByPos(i, j);
				if (config.script != "" && !config.isCloud) {
					auto zc = luaScript.runScript(config.script, z, c, config.maxIteration);
					z = std::get<0>(zc);
					c = std::get<1>(zc);
					n = std::get<2>(zc);
				}
				else {
					for (int k = 0; k < config.maxIteration; k++) {
						if (r > 0) {
							if (mod(z) > r) {
								kid = k;
								break;
							}
						}

						if (config.script == "") {
							z = _calcMondelbrot(z, config.dim) + c;
						}
						else {
							auto zc = luaScript.runScript(config.script, z, c, config.maxIteration);
							z = std::get<0>(zc);
							c = std::get<1>(zc);
						}
						if (config.isCloud) {
							int _x = z.real() * scaleX + (-(minX + range.x / 2) + config.width / 2);
							int _y = z.imag() * scaleY + (-(minY + range.y / 2) + config.height / 2);

							if (n < 10 || _x < 0 || _y < 0 || _x >= config.width || _y >= config.height) {
								n++;
								continue;
							}

							image.at<cv::Vec3b>(cv::Point(_x, _y)) =
								cv::Vec3b(claudColor.blue, claudColor.green, claudColor.red);
						}
					}
				}
				if (!config.isCloud) {
					double smoothed = log2(log2(z.real() * z.real() + z.imag() * z.imag()) / 2);
					int colorI = (int)(sqrt(kid - smoothed) * 256) % (int)(image.cols);
					//auto color = colorManager.getRgb(colorI);
					//auto color = defaultColorSheme(kid);
					auto color = colorManager.getColorByPos(colorI < 0 ? 0 :colorI, 0);//defaultColorSheme(kid);
					pixels[i][2] = color.red;
					pixels[i][1] = color.green;
					pixels[i][0] = color.blue;
				}
			}
		}
		threadsDone++;
	}


	DrawMachine dm;

	void RunLsystem() {
		initImage();
		auto startX = config.width / 2;
		auto startY = config.height / 2;
		dm.setStartXY(startX, startY);

		auto ls = config.getLsystem();
		ls.iterations = iterations;

		dm.runLsystem(ls);
		dm.draw(image, colorManager, config.outFilename, 
			config.isStepByStep, video, config.copyFrame);
	}

	void RunLapunov() {
		threadsDone = 0;
		maxX = centerX + range.x / 2;
		minX = centerX - range.x / 2;
		minY = centerY - range.y / 2;
		maxY = centerY + range.y / 2;
		dx = range.x / config.width;
		dy = range.y / config.height;

		int h = (int)(config.height / config.threadSize);
		int dh = config.height - h * config.threadSize;
		for (auto i = 0; i < config.threadSize; ++i) {
			std::thread t(&Fractal::RunLapunovThread, this, h * i, (i == config.threadSize - 1 ? h * (i + 1) + dh : h * (i + 1)));
			t.detach();
		}
		while (threadsDone != config.threadSize) {
			Sleep(10);
		}
	}


	void RunLapunovThread(int sj, int ej) {
		double x, y;
		register double xn, rn, lim=0.0;
		for (int j = sj; j < ej; j++) {
			y = minY + j * dy;
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i++) {
				x = maxX - i * dx;

				xn = config.xn;
				lim = 0.0;

				for (auto k = 0; k < config.maxIteration; k++) {
					auto id = k % config.lapunovPattern.size();
					if (config.lapunovPattern[id] == 'A') {
						rn = x;
					}
					else {
						rn = y;
					}
					xn = rn * xn * (1 - xn);
					lim += std::log(std::abs(rn * (1 - 2 * xn)));
				}

				lim /= 500;
				
				if (lim < -1000000000) {
					lim = -image.cols+1;
				}
				if (lim > 1000000000) {
					lim = image.cols-1;
				}
				if (lim >= config.maxIteration) {
					pixels[i][2] = 0;
					pixels[i][1] = 0;
					pixels[i][0] = 0;
				}
				else {
					int colorI = (int)(sqrt(lim*lim + 10) * 256) % image.cols;
					auto color = colorManager.getColorByPos(colorI < 0 ? 0 : colorI, 0);
					pixels[i][2] = color.red;
					pixels[i][1] = color.green;
					pixels[i][0] = color.blue;
				}
			}
		}
		threadsDone++;
	}

	std::vector<std::complex<double>> roots;

	std::complex<double> func(std::complex<double> in) {
		std::complex<double> res(0, 0);
		for (auto i = 0; i <= config.polynomDegree; i++) {
			res += config.coefs[i] * pow(in, config.polynomDegree - i);
		}
		return res;
	}

	std::complex<double> deriv(std::complex<double> in) {
		std::complex<double> res(0, 0);
		for (auto i = 0; i <= config.polynomDegree; i++) {
			if (config.polynomDegree - i == 0) {
				break;
			}
			res += config.coefs[i] * std::complex<double>(
				config.polynomDegree - i, 0) * pow(in, config.polynomDegree - i - 1);
		}
		return res;
	}
	
	void RunFractalNewton() {
		double err;

		std::vector<double> quad = {2.71828e-1, 3.14159e-1};
		std::vector<double> x;
		x.resize(config.polynomDegree + 1);
		GetQuads(config.coefs, config.polynomDegree, quad, x);
		int numr = GetRoots(x, config.polynomDegree, roots);

		threadsDone = 0;
		maxX = centerX + range.x / 2;
		minX = centerX - range.x / 2;
		minY = centerY - range.y / 2;
		maxY = centerY + range.y / 2;
		dx = range.x / config.width;
		dy = range.y / config.height;

		int h = (int)(config.height / config.threadSize);
		int dh = config.height - h * config.threadSize;
		for (auto i = 0; i < config.threadSize; ++i) {
			std::thread t(&Fractal::RunFractalNewtonTread, this, 
				h * i, (i == config.threadSize - 1 ? h * (i + 1) + dh : h * (i + 1)));
			t.detach();
		}
		while (threadsDone != config.threadSize) {
			Sleep(10);
		}
	}

	void RunFractalNewtonTread(int sj, int ej) {
		double x, y;
		register int n;
		double _c = cos(angle), _s = sin(angle);
		for (int j = sj; j < ej; j++) {
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i++) {
				x = maxX - i * dx;
				y = minY + j * dy;
				std::complex<double> z((x * _c - y * _s), (x * _s + y * _c));
				
				bool isSet = false;
				n = 0;
				while (n < config.maxIteration) {
					z -= func(z) / deriv(z);
					double tolerance = 0.000001;
					for (int r = 0; r < roots.size(); r++) {
						std::complex<double> difference = z - roots[r];
						if (fabs(std::real(difference)) < tolerance && fabs(std::imag(difference)) < tolerance) {
							isSet = true;

							//double smoothed = log2(log2(z.real() * z.real() + z.imag() * z.imag()) / 2);
							//int colorI = (int)(sqrt(n - smoothed) * 256) % (image.cols / config.polynomDegree);
							auto color = colorManager.getRgb((image.cols / config.polynomDegree)*r);

							pixels[i][2] = color.blue - (color.red != 0 ? (char)(n / (double)config.maxIteration * 2 * 128000) : 0);
							pixels[i][1] = color.green - (color.red != 0 ? (char)(n / (double)config.maxIteration * 2 * 128000) : 0);
							pixels[i][0] = color.red - (color.red != 0 ? (char)(n / (double)config.maxIteration * 2 * 128000) : 0);
							break;
						}
					}
					if (isSet) {
						break;
					}
					n++;
				}
				if (!isSet) {
					pixels[i][2] = 0;
					pixels[i][1] = 0;
					pixels[i][0] = 0;
				}
			}
		}
		threadsDone++;
	}

	void GetQuads(std::vector<double>& coefs,
		int polynomSize, std::vector<double>& quad,
		std::vector<double>& x) {
		double err, tmp;
		std::vector<double> b, z;
		double xr, xs;
		int iter, i, m;

		if ((tmp = coefs[0]) != 1.0) {
			coefs[0] = 1.0;
			for (i = 1; i <= polynomSize; i++) {
				coefs[i] /= tmp;
			}
		}
		if (polynomSize == 2) {
			x[0] = coefs[1];
			x[1] = coefs[2];
			return;
		}
		else if (polynomSize == 1) {
			x[0] = coefs[1];
			return;
		}
		m = polynomSize;
		b.resize(polynomSize + 1);
		z.resize(polynomSize + 1);
		b[0] = 1.0;
		for (i = 0; i <= polynomSize; i++) {
			z[i] = coefs[i];
			x[i] = 0.0;
		}
		do {
			if (polynomSize > m) {
				quad[0] = 3.14159e-1;
				quad[1] = 2.78127e-1;
			}
			do {// This loop tries to assure convergence
				for (i = 0; i < 5; i++) {
					FindQuad(z, m, b, quad, err, iter);
					if ((err > 1e-7) || (iter > config.maxIteration)) {
						DiffPoly(z, m, b);
						auto mm = m - 1;
						Recurse(z, m, b, mm, quad, err, iter);
					}
					Deflate(z, m, b, quad, err);
					if (err < 0.001) {
						break;
					}
					quad[0] = std::rand() % 8 - 4.0;
					quad[1] = std::rand() % 8 - 4.0;
				}
				if (err > 0.01) {
					std::cout << "Error! Convergence failure in quadratic x^2 + r*x + s." << std::endl;
					std::cout << "Enter new trial value for 'r': ";
					std::cin >> quad[1];
					std::cout << "Enter new trial value for 's' ( 0 to exit): ";
					std::cin >> quad[0];
					if (quad[0] == 0) {
						exit(1);
					}
				}
			} while (err > 0.01);
			x[m - 2] = quad[1];
			x[m - 1] = quad[0];
			m -= 2;
			for (i = 0; i <= m; i++) {
				z[i] = b[i];
			}
		} while (m > 2);
		if (m == 2) {
			x[0] = b[1];
			x[1] = b[2];
		}
		else {
			x[0] = b[1];
		}
	}

	int GetRoots(std::vector<double>& coefs, int& polynomSize,
		std::vector<std::complex<double>>& roots) {
		double sq, b2, c, disc;
		int i, m, numroots;
		roots.resize(polynomSize);
		m = polynomSize;
		numroots = 0;
		while (m > 1) {
			b2 = -0.5 * coefs[m - 2];
			c = coefs[m - 1];
			disc = b2 * b2 - c;
			if (disc < 0.0) {// complex roots
				sq = sqrt(-disc);
				roots[m - 2] = std::complex<double>(b2, sq);
				roots[m - 1] = std::complex<double>(b2, -sq);
				numroots += 2;
			}
			else {// real roots
				sq = sqrt(disc);
				double newRoot1 = fabs(b2) + sq;
				double newRoot2 = 0;
				if (b2 < 0.0) {
					newRoot1 = -newRoot1;
				}
				if (newRoot1 == 0) {
					newRoot2 = 0;
				}
				else {
					newRoot2 = c / newRoot1;
					numroots += 2;
				}
				roots[m - 2] = std::complex<double>(newRoot1, 0.0);
				roots[m - 1] = std::complex<double>(newRoot2, 0.0);
			}
			m -= 2;
		}
		if (m == 1) {
			roots[0] = std::complex<double>(-coefs[0], 0.0);
			numroots++;
		}
		return numroots;
	}

	void Deflate(std::vector<double>& coefs,
		int polynomSize,
		std::vector<double>& b, std::vector<double>& quad, double& err) {
		double r, s;
		int i;

		r = quad[1];
		s = quad[0];

		b[1] = coefs[1] - r;

		for (i = 2; i <= polynomSize; i++) {
			b[i] = coefs[i] - r * b[i - 1] - s * b[i - 2];
		}
		err = fabs(b[polynomSize]) + fabs(b[polynomSize - 1]);
	}

	void FindQuad(std::vector<double>& coefs,
		int polynomSize, std::vector<double>& b,
		std::vector<double>& quad,
		double& err, int& iter) {
		double dn, dr, ds, drn, dsn, eps, r, s;
		int i;

		std::vector<double> c;

		c.resize(polynomSize + 1);
		c[0] = 1.0;
		r = quad[1];
		s = quad[0];
		eps = 1e-15;
		iter = 1;

		do {
			if (iter > config.maxIteration) {
				break;
			}
			if (((iter) % 200) == 0) {
				eps *= 10.0;
			}
			b[1] = coefs[1] - r;
			c[1] = b[1] - r;

			for (i = 2; i <= polynomSize; i++) {
				b[i] = coefs[i] - r * b[i - 1] - s * b[i - 2];
				c[i] = b[i] - r * c[i - 1] - s * c[i - 2];
			}
			dn = c[polynomSize - 1] * c[polynomSize - 3] - c[polynomSize - 2] * c[polynomSize - 2];
			drn = b[polynomSize] * c[polynomSize - 3] - b[polynomSize - 1] * c[polynomSize - 2];
			dsn = b[polynomSize - 1] * c[polynomSize - 1] - b[polynomSize] * c[polynomSize - 2];

			if (fabs(dn) < 1e-10) {
				if (dn < 0.0) {
					dn = -1e-8;
				}
				else {
					dn = 1e-8;
				}
			}
			dr = drn / dn;
			ds = dsn / dn;
			r += dr;
			s += ds;
			iter++;
		} while ((fabs(dr) + fabs(ds)) > eps);
		quad[0] = s;
		quad[1] = r;
		err = fabs(ds) + fabs(dr);
	}

	void DiffPoly(std::vector<double>& coefs,
		int polynomSize, std::vector<double>& b) {
		double coef;
		int i;

		coef = (double)polynomSize;
		b[0] = 1.0;
		for (i = 1; i < polynomSize; i++) {
			b[i] = coefs[i] * ((double)(polynomSize - i)) / coef;
		}
	}

	void Recurse(std::vector<double>& coefs,
		int polynomSize, std::vector<double>& b, int& m,
		std::vector<double>& quad,
		double& err, int& iter) {
		double tst, e1, e2;

		std::vector<double> c, x, rs;
		rs.resize(2);

		if (fabs(b[m]) < 1e-16) {
			m--;
		}
		if (m == 2) {
			quad[0] = b[2];
			quad[1] = b[1];
			err = 0;
			iter = 0;
			return;
		}

		c.resize(m + 1);
		x.resize(polynomSize + 1);
		c[0] = x[0] = 1.0;
		rs[0] = quad[0];
		rs[1] = quad[1];
		iter = 0;
		FindQuad(b, m, c, rs, err, iter);
		tst = fabs(rs[0] - quad[0]) + fabs(rs[1] - quad[1]);
		if (err < 1e-12) {
			quad[0] = rs[0];
			quad[1] = rs[1];
		}
		if (((iter > 5) && (tst < 1e-4)) || ((iter > 20) && (tst < 1e-1))) {
			DiffPoly(b, m, c);
			int mm = m - 1;
			Recurse(coefs, polynomSize, c, mm, rs, err, iter);
			quad[0] = rs[0];
			quad[1] = rs[1];
		}
	}

	void Brownian() {
		int saveId = 1;
		auto startX = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / config.width));
		auto startY = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / config.height));
		for (int i = 0; i < config.iterations; i++) {
			auto rAngle = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * M_PI)));

			auto x = startX + config.distance * cos(rAngle);
			auto y = startY + config.distance * sin(rAngle);
			if (x < 0) x = 0;
			if (x >= config.width) x = config.width-1;
			if (y < 0) y = 0;
			if (y >= config.height) y = config.height - 1;

			auto _rgb = colorManager.getColorByPos(x, y);
			line(image,
				cv::Point(startX, startY),
				cv::Point(x, y),
				cv::Scalar(_rgb.blue, _rgb.green, _rgb.red),
				config.lineWidth,
				cv::LINE_AA);
			if (config.isStepByStep) {
				save(saveId); saveId++;
			}

			if (config.drawCircles) {
				auto _rgb2 = colorManager.getColorByPos(startX, startY);
				circle(image, cv::Point(x, y), config.drawCircles, cv::Scalar(_rgb.blue, _rgb.green, _rgb.red), cv::FILLED, 8, 0);
				circle(image, cv::Point(startX, startY), config.drawCircles, cv::Scalar(_rgb2.blue, _rgb2.green, _rgb2.red), cv::FILLED, 8, 0);
			}

			startX = x;
			startY = y;
		}
	}

	void DLA() {
		int saveId = 1;
		std::vector<std::pair<cv::Point, int>> circles;
		int maxX, minX, maxY, minY;

		auto checkX = [this](int& x, int _min, int _max) {
			if (x < _min) x = _min;
			if (x >= _max) x = _max - 1;
		};
		auto checkY = [this](int& y, int _min, int _max) {
			if (y < _min) y = _min;
			if (y >= _max) y = _max - 1;
		};

		auto circleCircleCollide = [this](int x, int y, int radius, std::pair<cv::Point, int> c) -> bool {
			return (c.first.x - x) * (c.first.x - x) + (y - c.first.y) * (y - c.first.y) <= (radius + c.second) * (radius + c.second);
		};
		
		auto startX = config.width / 2;
		auto startY =  config.height / 2;
		minX = startX - 10; checkX(minX, 0, config.width);
		maxX = startX + 10; checkX(maxX, 0, config.width);
		minY = startY - 10; checkY(minY, 0, config.height);
		maxY = startY + 10; checkY(maxY, 0, config.height);

		circles.push_back(std::make_pair(cv::Point(startX, startY), config.radius));

		auto _rgb = colorManager.getColorByPos(startX, startY);
		circle(image, cv::Point(startX, startY), config.radius, cv::Scalar(_rgb.blue, _rgb.green, _rgb.red), cv::FILLED, 8, 0);

		
		for (int i = 0; i < config.iterations; i++) {
			auto x = 0;//minX + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxX-minX-1)));
			auto y = 0;//minY + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxY-minY-1)));
			auto r = rand() % 4;
			if (r == 0) {
				x = minX;
				y = minY + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxY - minY - 1)));
			}
			else if (r == 1) {
				x = maxX;
				y = minY + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxY - minY - 1)));
			}
			else if (r == 2) {
				x = minX + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxX - minX - 1)));
				y = minY;
			}
			else if (r == 3) {
				x = minX + static_cast <int> (rand()) / (static_cast <int> (RAND_MAX / (maxX - minX - 1)));
				y = maxY;
			}
			while (true) {
				bool check = false;
				for (auto c : circles) {
					if (circleCircleCollide(x, y, config.radius, c)) {
						check = true;
						break;
					}
				}
				if (check) {
					std::cout << i << std::endl;
					checkX(x, 0, config.width-1);
					checkY(y, 0, config.height-1);
					auto _rgb = colorManager.getColorByPos(x, y);
					circle(image, cv::Point(x, y), config.radius, cv::Scalar(_rgb.blue, _rgb.green, _rgb.red), cv::FILLED, 8, 0);
					if (config.isStepByStep) {
						save(saveId); saveId++;
					}
					if (x + 5 > maxX) maxX = x + 5;
					if (x - 5 < minX) minX = x - 5;
					if (y + 5 > maxY) maxY = y + 5;
					if (y - 5 < minY) minY = y - 5;
					checkX(maxX, 0, config.width);
					checkY(maxY, 0, config.height);
					checkX(minX, 0, config.width);
					checkY(minY, 0, config.height);
					circles.push_back(std::make_pair(cv::Point(x, y), config.radius));
					break;
				}
				x = x + -1 + rand() % 3;
				y = y + -1 + rand() % 3;
				checkX(x, minX, maxX);
				checkY(y, minY, maxY);

			}
		}

	}

	void Attractor() {
		LuaScript luaScript;
		if (config.script != "") {
			luaScript.init();
		}
		else {
			E_LOG("Attractor: Path to script empy.");
			return;
		}
		luaScript.runAttractorScript(config.script, config.iterations, config.width, config.height, config.scale);
		auto num = 0;
		std::ifstream fin("attractorGrid.txt");
		for (int j = 0; j < config.height; j ++) {
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i++) {
				fin >> num;
				if (num == 0) {
					continue;
				}
				num = image.cols * num / config.colorScale;
				int colorI = (int)(sqrt(num + 10) * 256) % image.cols;
				auto color = colorManager.getColorByPos(colorI < 0 ? 0 : colorI, 0);
				pixels[i][2] = color.red;
				pixels[i][1] = color.green;
				pixels[i][0] = color.blue;
			}
			fin >> num;
		}
	}


	int recursion = 0;
	int pen = 0;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution;

	void Lightning() {
		recursion = 0;
		auto startY = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / config.height / 2));
		auto seed = config.seed;
		if (seed < 0) seed = std::chrono::system_clock::now().time_since_epoch().count();
		generator = std::default_random_engine(seed);
		lightningFun(std::vector<cv::Point>{cv::Point(0, config.height / 4 + startY)}, std::vector<double>{(5.0 + rand()%30)}, config.spred);
		dm.draw(image, colorManager, config.outFilename,
			config.isStepByStep, video, config.copyFrame);
	}

	double clamp(double n, double minn, double maxn) {
		return max(min(maxn, n), minn);
	}

	void lightningFun(std::vector<cv::Point> pos, std::vector<double> head, int spread) {
		auto x1 = 0;
		auto y1 = 0;
		auto h3 = 0;

		auto x2 = 0;
		auto y2 = 0;
		auto h4 = 0;
		int saveId = 1;
		recursion++;
		auto ps = pos.size();
		for (auto i = 0; i < ps; i++) {
			auto _pen = config.penSize - clamp(2 * std::floor(std::log2(pos.size())), 0, config.penSize-1);

			dm.setStartXY(pos[i].x, pos[i].y);

			for (auto j = 0; j < 20 + rand() % 20; j++) {
				pen = _pen;
				distribution = std::normal_distribution<double>(spread, 15);
				auto head1 = head[i] + distribution(generator);
				if (180 > head1 > 0) {
					dm.setAlpha(-head1);
				}
				else {
					dm.setAlpha(head1);
				}
				distribution = std::normal_distribution<double>(1, 4);
				dm.forward(pen*std::abs(std::round(distribution(generator))));
				if (config.isStepByStep) {
					save(saveId); saveId++;
				}
				x1 = dm.getX();
				y1 = dm.getY();
				h3 = dm.getA();
				
			}

			dm.setStartXY(pos[i].x, pos[i].y);
			pos.push_back(cv::Point(x1, y1));
			head.push_back(h3);

			for (auto k = 0; k < 20 + rand() % 20; k++) {
				pen = _pen;

				distribution = std::normal_distribution<double>(spread, 15);
				auto head2 = (head[i] + distribution(generator));
				if (180 > head2 > 0) {
					dm.setAlpha(-head2);
				}
				else {
					dm.setAlpha(head2);
				}

				distribution = std::normal_distribution<double>(2, 4);
				dm.forward(pen * abs(round(distribution(generator))));
				if (config.isStepByStep) {
					save(saveId); saveId++;
				}
				x2 = dm.getX();
				y2 = dm.getY();
				h4 = dm.getA();
				
			}
			pos.push_back(cv::Point(x2, y2));
			head.push_back(h4);

			auto _min = 100000000000;
			for (auto p : pos) {
				auto x2y2 = p.x * p.x + p.y * p.y;
				if (x2y2 < _min) _min = x2y2;
			}
			
			if (recursion < 18 && std::sqrt(_min) < 50000) {
				lightningFun(std::vector<cv::Point>{pos.begin() + (i + 1), pos.end()}, std::vector<double>{head.begin() + (i + 1), head.end()}, spread);
			}
		}
	}

	std::vector<std::vector<double>> terrainArray;
	
	void generateTerrain(int leftIndex, int rightIndex, double displacement) {

		if ((leftIndex + 1) == rightIndex) return;
		auto midIndex = floor((leftIndex + rightIndex) / 2);
		auto change = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 2 - 1) * displacement;
		terrainArray[0][midIndex] = (terrainArray[0][leftIndex] + terrainArray[0][rightIndex]) / 2 + change;

		displacement = displacement * config.roughness;
		generateTerrain(leftIndex, midIndex, displacement);
		generateTerrain(midIndex, rightIndex, displacement);
	}

	
	void Landscape2d() {
		terrainArray = std::vector<std::vector<double>>(config.height, std::vector<double>(config.width, 0));

		generateTerrain(0, config.width-1, config.displacement);

		for (auto i = 0; i < terrainArray[0].size(); ++i) {
			auto _rgb = colorManager.getColorByPos(i, 0);
			line(image,
				cv::Point(i, 0),
				cv::Point(i, (terrainArray[0][i]+1.0) * config.lsScale),
				cv::Scalar(_rgb.blue, _rgb.green, _rgb.red),
				config.lineWidth,
				cv::LINE_AA);
		}
		cv::Point2f srcCenter(image.cols / 2.0, image.rows / 2.0);
		auto rotMat = cv::getRotationMatrix2D(srcCenter, 180, 1.0);
		warpAffine(image, image, rotMat, image.size());
	}

	int levels = 10;
	int divisions = 0;

	double _min =0.0, _max = 0.0;
	void Landscape3d() {
		divisions = max(config.width, config.height); //1 << levels;
		levels = log2(divisions)+1;
		divisions = 1 << levels;
		terrainArray = std::vector<std::vector<double>>(divisions+1, std::vector<double>(divisions+1, 0));

		terrainArray[0][0] = rnd();
		terrainArray[0][divisions] = rnd();
		terrainArray[divisions][divisions] = rnd();
		terrainArray[divisions][0] = rnd();
		double rough = config.roughness;
		for (int i = 0; i < levels; ++i) {
			int q = 1 << i, r = 1 << (levels - i), s = r >> 1;
			for (int j = 0; j < divisions; j += r)
				for (int k = 0; k < divisions; k += r)
					diamond(j, k, r, rough);
			if (s > 0)
				for (int j = 0; j <= divisions; j += s)
					for (int k = (j + s) % r; k <= divisions; k += r)
						square(j - s, k - s, r, rough);
			rough *= config.roughness;
		}
		_min = _max = terrainArray[0][0];
		for (int i = 0; i <= divisions; ++i)
			for (int j = 0; j <= divisions; ++j)
				if (terrainArray[i][j] < _min) _min = terrainArray[i][j];
				else if (terrainArray[i][j] > _max) _max = terrainArray[i][j];

		for (int j = 0; j < config.height; j++) {
			cv::Vec3b* pixels = image.ptr<cv::Vec3b>(j);
			for (int i = 0; i < config.width; i++) {

				int colorI = (int)(sqrt(terrainArray[j][i] * config.lsScale + 10) * 256) % image.cols;
				auto color = colorManager.getColorByPos(colorI < 0 ? 0 : colorI, 0);
				pixels[i][2] = color.red;
				pixels[i][1] = color.green;
				pixels[i][0] = color.blue;
			}
		}
	}

	void diamond(int x, int y, int side, double scale) {
		if (side > 1) {
			int half = side / 2;
			double avg = (terrainArray[x][y] + terrainArray[x + side][y] +
				terrainArray[x + side][y + side] + terrainArray[x][y + side]) * 0.25;
			terrainArray[x + half][y + half] = avg + rnd() * scale;
		}
	}
	void square(int x, int y, int side, double scale) {
		int half = side / 2;
		double avg = 0.0, sum = 0.0;
		if (x >= 0) {
			avg += terrainArray[x][y + half]; sum += 1.0;
		}
		if (y >= 0) {
			avg += terrainArray[x + half][y]; sum += 1.0;
		}
		if (x + side <= divisions) {
			avg += terrainArray[x + side][y + half]; sum += 1.0;
		}
		if (y + side <= divisions) {
			avg += terrainArray[x + half][y + side]; sum += 1.0;
		}
		terrainArray[x + half][y + half] = avg / sum + rnd() * scale;
	}
	double rnd() {
		return 2.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX) - 1.0;
	}
	double getAltitude(double i, double j) {
		double alt = terrainArray[(int)(i * divisions)][(int)(j * divisions)];
		return (alt - _min) / (_max - _min);
	}

	double getProbSum() {
		return 1;
	}

	
	void Fern() {
		auto checkX = [this](int x, int _min, int _max) {
			if (x < _min) return _min;
			if (x >= _max) return _max - 1;
			return x;
		};
		auto checkY = [this](int y, int _min, int _max) {
			if (y < _min) return  _min;
			if (y >= _max) return _max - 1;
			return y;
		};

		int saveId = 1;
		
		double x = 0;
		double y = 0;
		for (auto n = 0; n < config.iterations; n++) {
			auto _x = checkX(x*config.coordScaleX + config.width / 2, 0, config.width);
			auto _y = checkY(y*config.coordScaleY, 0, config.height);
			auto _rgb = colorManager.getColorByPos(_x, _y);
			circle(image, cv::Point(_x, _y), 3, cv::Scalar(_rgb.blue, _rgb.green, _rgb.red), cv::FILLED, 8, 0);
			if (config.isStepByStep) {
				save(saveId); saveId++;
			}
			r = rand() / static_cast<double>(RAND_MAX)  * getProbSum();
			if (r == 0) {
				r += 0.0000001;
			}
			double xn = x;
			double yn = y;

			auto val = 0.0;
			int id = 0;
			while (r > val) {
				val += config.fs[id][config.fs[0].size() - 1];
				id++;
			}
			id--;
			x = xn * config.fs[id][0] + yn * config.fs[id][1] + config.fs[id][4];
			y = xn * config.fs[id][2] + yn * config.fs[id][3] + config.fs[id][5];
		}

		cv::Point2f srcCenter(image.cols / 2.0, image.rows / 2.0);
		cv::Mat rotMat = getRotationMatrix2D(srcCenter, 180, 1.0);
		warpAffine(image, image, rotMat, image.size());
	}
};
