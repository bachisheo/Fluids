#include "ImguiWindow.h"

#include <SFML/Window/Event.hpp>

#include "imgui-SFML.h"
#include "imgui.h"

ImguiWindow::ImguiWindow(FluidCube* fluid)
{
	_window.create(sf::VideoMode(400, 400), "");
	_window.setPosition({0,0});
	if(!ImGui::SFML::Init(_window))
		throw std::exception("Окно ImGui не загрузилось");
	_fluid = fluid;
	_diff = fluid->GetDiffusionCoef();
	_visc = fluid->GetViscosityCoef();

}

void ImguiWindow::update(float dt)
{
	ImGui::SFML::Update(_window, sf::seconds(dt));

	ImGui::Begin("Editor");
	if (ImGui::DragFloat("Diffusion", &_diff, 0.1f, 0, 100))
		_fluid->SetDiffusionCoef(_diff);
	if (ImGui::DragFloat("Viscosity", &_visc, 0.1f, 0, 100))
		_fluid->SetViscosityCoef(_visc);
	if (ImGui::DragInt("Cursor Radius", &_radius, 0.1, 0, 20))
		_fluid->SetCursorRadiusCoef(_radius);
	ImGui::End();
}

void ImguiWindow::pollEvents()
{
	sf::Event event{};
	while (_window.pollEvent(event))
		ImGui::SFML::ProcessEvent(_window, event);
}

void ImguiWindow::draw()
{
	_window.clear();
	ImGui::SFML::Render(_window);
	_window.display();
}
