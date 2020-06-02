#pragma once
#include <cstdlib>
#include <ctime>
#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>

struct rgb {
    rgb(int r, int g=0, int b=0) : red(r), green(g), blue(b) {}
    rgb(std::vector<int> vec) {
        setColorVec(vec);
    }
    rgb() {}
    int red = 0;
    int green = 0;
    int blue = 0;

    bool operator==(rgb lhr) {
        return red == lhr.red && green == lhr.green && blue == lhr.blue;
    }

    std::vector<int> getColorVec() {
        return std::vector<int> {red, green, blue};
    }

    void setColorVec(std::vector<int> vec) {
    	if (vec.size() != 3) {
    		//ERROR
            return;
    	}
        red = vec[0];
        green = vec[1];
        blue = vec[2];
    }
	
    static rgb genRand() {
        return rgb(rand() % 255, rand() % 255, rand() % 255);
    }
};

static rgb INVALID_COLOR = {0, 0, 0};

class ColorManager {
public:
    ColorManager(int seed=-1) {
        if (seed < 0)
            srand(time(NULL));
        else
            srand(seed);
    }

    ColorManager(int min, int max, std::vector<rgb> stops,
        rgb minOutlierColor = INVALID_COLOR, rgb maxOutlierColor = INVALID_COLOR, int seed = -1) {
        if (seed < 0)
            srand(time(NULL));
        else
            srand(seed);
    	initialize(min, max, stops, minOutlierColor, maxOutlierColor);
    }

    virtual ~ColorManager() {}

    void initialize(int min,
        int max,
        std::vector<rgb> stops,
        rgb minOutlierColor = INVALID_COLOR,
        rgb maxOutlierColor = INVALID_COLOR) {
        _min = min;
        _max = max;
        _stops = stops;
        _minOutlierColor = minOutlierColor;
        _maxOutlierColor = maxOutlierColor;
    }

    rgb getRgb(int value) {
        if (value < _min) {
            return _minOutlierColor ==
                INVALID_COLOR ? _stops.front() : _minOutlierColor;
        }
        if (value > _max) {
            return _maxOutlierColor ==
                INVALID_COLOR ? _stops.back() : _maxOutlierColor;
        }

        int range = _max - _min;
        int v = value - _min;
        float step = range / (float)(_stops.size() - 1);
        int bin = (int)(v / step);

        float normalized_v = (v - bin * step) / step;

        return interpolate(_stops[bin], _stops[bin + 1], normalized_v);
    }


    cv::Mat gradImg;
    void genColorGradientMat(std::vector<std::vector<rgb>>& colorsMat) {
        auto stepX = 60;
        auto stepSizeX = 1;
        auto stepY = 50; //height
        auto stepSizeY = 1; //height

        std::vector<std::vector<std::vector<int>>> rows;
        for (auto j = 0; j < colorsMat.size(); j++) {
            std::vector<std::vector<int>> row;
        	if (colorsMat[j].size() == 1) {
                colorsMat[j].push_back(colorsMat[0][0]);
        	}
            for (auto i = 0; i < colorsMat[j].size() - 1; i++) {
                auto vec = gradient(colorsMat[j][i], colorsMat[j][i + 1], stepX, stepSizeX);
                row.insert(row.end(), vec.begin(), vec.end());
            }
            if (row.size() == 1) {
                row.push_back(row[0]);
            }
            rows.push_back(row);
        }
        if (rows.size() == 1) {
            rows.push_back(rows[0]);
        }

        if (rows.size() == 0) {
            return;
        }
    	
        gradImg = cv::Mat::zeros(stepY * stepSizeY * (rows.size() - 1), rows[0].size(), CV_8UC3);
        cv::rotate(gradImg, gradImg, cv::ROTATE_90_CLOCKWISE);

        for (auto k = 0; k < rows.size() - 1; k++) {
            auto top_color = rows[k];
            auto bottom_color = rows[k + 1];
            auto start = k * stepY * stepSizeY;
            auto end = (k + 1) * stepY * stepSizeY;
            for (auto i = 0; i < rows[0].size(); i++) {
                cv::Vec3b* pixels = gradImg.ptr<cv::Vec3b>(i);
                auto col = gradient(rgb(top_color[i][0], top_color[i][1], top_color[i][2]),
                    rgb(bottom_color[i][0], bottom_color[i][1], bottom_color[i][2]), stepY, stepSizeY);
                for (auto j = start; j < end; j++) {
                    pixels[j][0] = col[j - start][2];
                    pixels[j][1] = col[j - start][1];
                    pixels[j][2] = col[j - start][0];
                }
            }
        }

        cv::rotate(gradImg, gradImg, cv::ROTATE_90_COUNTERCLOCKWISE);

        initialize(0, 700, colorsMat[0]);
    }

	//Use image to get color
    void genColorMatFromImage(std::string& path) {
        gradImg = cv::imread(path, cv::COLOR_RGB2BGR);
    }

	void resizeColorGrad(int w, int h) {
    	if (gradImg.cols == 0 && gradImg.rows == 0) {
            gradImg = cv::Mat::zeros(h, w, CV_8UC3);
    		return;
    	}
        cv::resize(gradImg, gradImg, cv::Size(w, h), 0.0, 0.0);
    }

	rgb getColorByPos(int x, int y) {
    	auto c = gradImg.at<cv::Vec3b>(y, x);
        return rgb(c[2], c[1], c[0]);
    }

    std::vector<std::vector<int>> gradient(rgb start, rgb end, float steps, int stepSize) {
        float diff_r = ((end.red - start.red) / steps);
        float diff_g = ((end.green - start.green) / steps);
        float diff_b = ((end.blue - start.blue) / steps);

        std::vector<std::vector<int>> colors;
        for (auto k = 0; k < stepSize; k++)
            colors.push_back(std::vector<int>{start.red, start.green, start.blue});
        for (auto i = 0; i < steps - 2; i++) {
            std::vector<int> color;
            color.push_back(static_cast<int>(start.red + i * diff_r));
            color.push_back(static_cast<int>(start.green + i * diff_g));
            color.push_back(static_cast<int>(start.blue + i * diff_b));
            for (auto k = 0; k < stepSize; k++) {
                colors.push_back(color);
            }
        }
        for (auto k = 0; k < stepSize; k++)
            colors.push_back(std::vector<int>{end.red, end.green, end.blue});
        return colors;
    }
private:
    int _min = 0;
    rgb _minColor;
    rgb _minOutlierColor;

    int _max = 0;
    rgb _maxColor;
    rgb _maxOutlierColor;

    std::vector<rgb> _stops;

    rgb interpolate(rgb c1, rgb c2, float normalized_value) {
        if (normalized_value <= 0.0) { return c1; }
        if (normalized_value >= 1.0) { return c2; }

        int red = (int)((1.0 - normalized_value) * c1.red +
            normalized_value * c2.red);
        int green = (int)((1.0 - normalized_value) * c1.green +
            normalized_value * c2.green);
        int blue = (int)((1.0 - normalized_value) * c1.blue +
            normalized_value * c2.blue);

        return rgb(red, green, blue);
    };
};
