#include "3rd/argParse.h"
#include "fractal.h"
#include "config.h"

int main(int argc, char* argv []) {
	std::string configPath = "configs/config.toml";
	ap::ArgParse ap(L"MathOfArt");
	ap.addArg(std::vector<std::wstring>{L"-c", L"--config"}, ap::PayloadType::TYPE_WSTRING, [&ap, &configPath]() {
		auto res = ap.payloadMap.find(L"--config");
		configPath = ws2s((*(*res).second).wstrVal);
		std::wcout << "Config: " << (*(*res).second).wstrVal << std::endl;;
	}, new ap::Payload(""), L"config path");
	ap.parse(argc, argv);

	if (configPath == "") {
        E_LOG("Need config file");
        return -1;
	}
	
	Fractal f(configPath);
	f.run();
	//f.show();
	return 0;
}