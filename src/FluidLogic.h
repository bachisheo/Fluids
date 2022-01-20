#pragma once
#include <SFML/Graphics/Image.hpp>
#include "Field.h"
class FluidLogic
{
public:
	virtual ~FluidLogic() = default;
	//the public functions
	virtual void AddSourceOfDensity(int x0, int y0, int x, int y) = 0;
	virtual void update() = 0;
	virtual void reset() = 0;
	virtual void render(int size, sf::Image& image) = 0;

protected:
	virtual void diffuseVelocity() = 0;
	virtual void diffuseDensity() = 0;
	virtual void advect() = 0;
	virtual void projectVelocity() = 0;
	virtual void setBoundary(int d, Field x) = 0;

	//void render(int sixe, sf::Image& image);
};
