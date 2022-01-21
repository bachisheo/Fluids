#include <functional>
#include <iostream>

#include "Fluid.h"
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/Image.hpp>
#include <SFML/System/Clock.hpp>

#include "DummyFluid.h"

void drawLine(sf::Vector2i p0, sf::Vector2i p1, std::function<void(int x, int y, int dx, int dy)> handler)
{
	auto p = p0;
	auto fdx = p1.x - p0.x;
	auto fdy = p1.y - p0.y;
	auto dx = abs(p0.x - p1.x);
	auto sx = p0.x < p1.x ? 1 : -1;
	auto dy = -abs(p1.y - p0.y);
	auto sy = p0.y < p1.y ? 1 : -1;
	auto err = dx + dy;
	while (true)
	{
		handler(p.x, p.y, fdx, fdy);
		if (p == p1) break;
		auto e2 = 2 * err;
		if (e2 >= dy)
		{
			err += dy;
			p.x += sx;
		}
		if (e2 <= dx)
		{
			err += dx;
			p.y += sy;
		}
	}
}
int main()
{
	//fluid parameters
	int xSize, ySize;
	xSize = ySize = 300;
	int elementSize = 3;
	float diff = 0.05;
	float visc = 0.4975;
	float dens_value = 20;
	//simulation parameters
	int numIterations = 5;
	float eps = 1;
	//mouse parameters
	int mouseX = 0;
	int mouseY = 0;

	//setting up render window
	sf::RenderWindow window(sf::VideoMode((xSize)*elementSize + 1, (ySize)*elementSize + 1), "Project");


	//creating rendering tools
	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;
	image.create((xSize)*elementSize + 1, (ySize)*elementSize + 1, sf::Color::Black);

	//now let's create the fluid we want (x, y, diffusion, viscosity, timestep, accuracy/quality)
	auto fluid = new FluidCube(ySize, diff, visc,  numIterations);
	sf::Vector2i prevPosition = sf::Mouse::getPosition(window);
	int lastMouseX = prevPosition.x;
	int lastMouseY = prevPosition.y;
	bool is_pressed = 0;
	sf::Clock clock;
	while (window.isOpen())
	{
		auto dt = clock.getElapsedTime().asSeconds();
		if (dt == 0) continue;
		clock.restart();

		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
			if (event.type == sf::Event::MouseMoved) {
				lastMouseX = mouseX;
				lastMouseY = mouseY;
				sf::Vector2i localPosition = { event.mouseMove.x, event.mouseMove.y };
				prevPosition = { lastMouseX,lastMouseY };
				mouseX = localPosition.x;
				mouseY = localPosition.y;
		
				drawLine(prevPosition, localPosition, [&fluid, elementSize, dt](int x, int y, int dx, int dy) { fluid->addVelocity(x / elementSize, y / elementSize, dx / elementSize, dy / elementSize, dt); });

			}
			if(event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left)
			{
				is_pressed = true;
			}
			if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left)
			{
				is_pressed = false;
			}
			
		}
		if (is_pressed)
		{

			fluid->addDensity(lastMouseX / elementSize, lastMouseY / elementSize, dens_value);
			drawLine({mouseX, mouseY}, { mouseX, mouseY-10 }, [&fluid, elementSize, dt](int x, int y, int dx, int dy) { fluid->addVelocity(x / elementSize, y / elementSize, dx / elementSize, dy / elementSize, dt); });
			drawLine({ lastMouseX, lastMouseY }, { mouseX, mouseY },
				[&fluid, elementSize, lastMouseX, lastMouseY, dens_value](int x, int y, int dx, int dy) { fluid->addDensity(x / elementSize, y / elementSize, dens_value); });
		}

		
		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
			fluid->reset();
		}
		window.clear();
		fluid->update(dt);
		fluid->render(elementSize, image);
		texture.loadFromImage(image);
		sprite.setTexture(texture, true);
		window.draw(sprite);
		window.display();
	}
	return 0;
}
