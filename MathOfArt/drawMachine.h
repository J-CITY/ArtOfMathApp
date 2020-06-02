#pragma once

#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <stack>
#include <tuple>
#include "lsystem.h"
#include "utils.h"

struct DrawMachine {
private:
	int x = 0, y = 0;
	int alpha = 270;
	int polyId = -1;
	
	std::vector<std::tuple<cv::Point, cv::Point, int>> lines;
	std::vector<std::pair<cv::Point, int>> points;
	std::stack<std::tuple<int, int, int, float>> drawState;
	std::vector<std::vector<cv::Point>> polygons;
public:
	void setStartXY(int _x, int _y) {
		x = _x;
		y = _y;
	}
	
	void right(int in) {
		alpha += in;
		alpha = alpha % 360;
	}

	void left(int in) {
		alpha -= (in % 360);
		if (alpha < 0)
			alpha = 360 - (-1 * alpha);
	}

	void setAlpha(int in) {
		alpha = in;
	}

	int getX() { return x; }
	int getY() { return y; }
	int getA() { return alpha; }
	
	void forward(int sz) {
		auto xe = round((double)x + std::cos((double)alpha * M_PI / 180.0) * (double)sz);
		auto ye = round((double)y + std::sin((double)alpha * M_PI / 180.0) * (double)sz);

		//line(image,
		//     cv::Point(x, y),
		//     cv::Point(xe, ye),
		//     cv::Scalar(255, 0, 0),
		//	1,
		//     cv::LINE_AA);
		if (polyId < 0) {
			lines.push_back(std::tuple<cv::Point, cv::Point, int>(cv::Point(x, y), cv::Point(xe, ye), _ls.widthLine));
		}
		else {
			polygons[polyId].push_back(cv::Point(xe, ye));
		}
		x = xe;
		y = ye;
	}
	void forwardWithoutDraw(int sz) {
		auto xe = round((double)x + std::cos((double)alpha * M_PI / 180.0) * (double)sz);
		auto ye = round((double)y + std::sin((double)alpha * M_PI / 180.0) * (double)sz);
		x = xe;
		y = ye;
	}
	Lsystem _ls;
	void runLsystem(Lsystem& ls) {

		lines.clear();
		points.clear();
		drawState = std::stack<std::tuple<int, int, int, float>>();
		polygons.clear();
		_ls = ls;
		auto instructions = ls.calcLsystem();
		for (auto& ch : instructions) {
			switch (ch) {
			case 'F':
				forward(ls.distance);
				break;
			case '+':
				if (ls.swapPlusMinus) {
					right(ls.angle);
				}
				else {
					left(ls.angle);
				}
				break;
			case '-':
				if (ls.swapPlusMinus) {
					left(ls.angle);
				}
				else {
					right(ls.angle);
				}
				break;
			case 'f':
				forward(ls.distance);
				break;
			case '|':
				left(180);
				break;
			case '[':
				drawState.push(std::tuple<int, int, int, float>(x, y, alpha, ls.distance));
				break;
			case ']': {
				auto pos = drawState.top();
				drawState.pop();
				x = std::get<0>(pos);
				y = std::get<1>(pos);
				alpha = std::get<2>(pos);
				ls.distance = std::get<3>(pos);
				break;
			}
			case '>':
				ls.distance *= ls.scaleFactor;
				break;
			case '<':
				ls.distance /= ls.scaleFactor;
				break;
			case '#':
				ls.widthLine += ls.widthStep;
				break;
			case '!':
				ls.widthLine -= ls.widthStep;
				if (ls.widthLine < 1) {
					ls.widthLine = 1;
				}
				break;
			case '@':
				points.push_back(std::make_pair(cv::Point(x, y), ls.widthLine));
				break;
			case '{':
				polyId++;
				polygons.push_back(std::vector<cv::Point>());
				polygons[polyId].push_back(cv::Point(x, y));
				break;
			case '}':
				polyId--;
				break;
			case '&':
				ls.swapPlusMinus = !ls.swapPlusMinus;
				break;
			case '(':
				ls.angle += ls.angleStep;
				break;
			case ')':
				ls.angle -= ls.angleStep;
				break;
			default:
				break;
			}
		}
	}

	cv::Point shift;

	bool checkPoint(cv::Point& p, int width, int height) {
		if (p.x < 0 || p.x >= width || p.y < 0 || p.y >= height) {
			return false;
		}
		return true;
	}
	
	void draw(cv::Mat& img, ColorManager& cm, std::string saveName="", 
		bool needSave=false, cv::VideoWriter* video=nullptr, int copyFrame=10) {
		//calculatePointShift(img.cols, img.rows);
		int id = 0;
		for (auto& l : lines) {
			std::get<0>(l).x += shift.x;
			std::get<0>(l).y += shift.y;
			std::get<1>(l).x += shift.x;
			std::get<1>(l).y += shift.y;
			if (checkPoint(std::get<0>(l), img.cols, img.rows) && checkPoint(std::get<1>(l), img.cols, img.rows)) {
				auto _rgb = cm.getColorByPos(std::get<0>(l).x + (std::get<1>(l).x - std::get<0>(l).x) / 2,
					std::get<0>(l).y + (std::get<1>(l).y - std::get<0>(l).y) / 2);
				line(img,
					std::get<0>(l),
					std::get<1>(l),
					cv::Scalar(_rgb.blue, _rgb.green, _rgb.red),
					std::get<2>(l),
					cv::LINE_AA);
			}
			else if (checkPoint(std::get<0>(l), img.cols, img.rows) && !checkPoint(std::get<1>(l), img.cols, img.rows)) {
				auto _rgb = cm.getColorByPos(std::get<0>(l).x, std::get<0>(l).y);
				line(img,
					std::get<0>(l),
					std::get<1>(l),
					cv::Scalar(_rgb.blue, _rgb.green, _rgb.red),
					std::get<2>(l),
					cv::LINE_AA);
			}
			else if (!checkPoint(std::get<0>(l), img.cols, img.rows) && checkPoint(std::get<1>(l), img.cols, img.rows)) {
				auto _rgb = cm.getColorByPos(std::get<1>(l).x, std::get<1>(l).y);
				line(img,
					std::get<0>(l),
					std::get<1>(l),
					cv::Scalar(_rgb.blue, _rgb.green, _rgb.red),
					std::get<2>(l),
					cv::LINE_AA);
			}

			if (needSave) {
				if (video) {
					for (auto i = 0; i < copyFrame; i++) {
						video->write(img);
					}
				}
				else {
					cv::imwrite(saveName + std::to_string(id) + ".png", img);
				}
				id++;
			}
		}

		id = 0;
		for (auto& l : points) {
			l.first.x += shift.x;
			l.first.y += shift.y;
			auto _rgb = cm.getColorByPos(l.first.x, l.first.y);
			circle(img, cv::Point(l.first.x, l.first.y), l.second, cv::Scalar(_rgb.blue, _rgb.green, _rgb.red), cv::FILLED, 8, 0);
			if (needSave) {
				if (video) {
					for (auto i = 0; i < copyFrame; i++) {
						video->write(img);
					}
				}
				else {
					cv::imwrite(saveName + std::to_string(id) + ".png", img);
				}
				id++;
			}
		}
		id = 0;
		for (auto& l : polygons) {
			for (auto& p : l) {
				p.x += shift.x;
				p.y += shift.y;
			}
			rgb pcolor;
			for (auto& p : l) {
				if (checkPoint(p, img.cols, img.rows)) {
					pcolor = cm.getColorByPos(p.x, p.y);
					break;
				}
			}
			const cv::Point* elementPoints[1] = {&l[0]};
			int numberOfPoints = (int)l.size();
			cv::fillPoly(img, elementPoints, &numberOfPoints, 1, cv::Scalar(), 8);
			if (needSave) {
				if (video) {
					for (auto i = 0; i < copyFrame; i++) {
						video->write(img);
					}
				}
				else {
					cv::imwrite(saveName + std::to_string(id) + ".png", img);
				}
				id++;
			}
		}
	}

	void calculatePointShift(int width, int height) {
		shift.x = 0;
		shift.y = 0;

		int minX = 10000, maxX=0, minY = 10000, maxY=0;

		for (auto& l : lines) {
			if (std::get<0>(l).x < minX) minX = std::get<0>(l).x;
			if (std::get<0>(l).x > maxX) maxX = std::get<0>(l).x;
			if (std::get<1>(l).x < minX) minX = std::get<1>(l).x;
			if (std::get<1>(l).x > maxX) maxX = std::get<1>(l).x;
			
			if (std::get<0>(l).y < minY) minY = std::get<0>(l).y;
			if (std::get<0>(l).y > maxY) maxY = std::get<0>(l).y;
			if (std::get<1>(l).y < minY) minY = std::get<1>(l).y;
			if (std::get<1>(l).y > maxY) maxY = std::get<1>(l).y;
		}

		for (auto& l : points) {
			if (l.first.x < minX) minX = l.first.x;
			if (l.first.x > maxX) maxX = l.first.x;

			if (l.first.y < minY) minY = l.first.y;
			if (l.first.y > maxY) maxY = l.first.y;
		}

		for (auto& l : polygons) {
			for (auto& p : l) {
				if (p.x < minX) minX = p.x;
				if (p.x > maxX) maxX = p.x;

				if (p.y < minY) minY = p.y;
				if (p.y > maxY) maxY = p.y;
			}
		}

		auto cx = width / 2;
		auto cy = height / 2;

		auto cx2 = minX + (maxX - minX) / 2;
		auto cy2 = minY + (maxY - minY) / 2;

		shift.x = -cx2 + cx;
		shift.y = -cy2 + cy;
	}
};
