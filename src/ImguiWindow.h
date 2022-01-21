#pragma once
#include <SFML/Graphics/RenderWindow.hpp>
#include "DummyFluid.h"

class ImguiWindow 
{
public:
	explicit ImguiWindow(FluidCube* fluid);
	void update(float dt);
	void setPosition(const sf::Vector2i& position) { this->_position = position; }
	void pollEvents();
	void draw();

private:

	sf::RenderWindow _window;
	float _diff, _visc;
	int _radius;
	sf::Vector2i _position;
	FluidCube* _fluid;
};
