#pragma once

#include <UGM/UGM.h>

struct CanvasData {
	std::vector<Ubpa::pointf2> points;
	std::vector<Ubpa::pointf2> curve1;
	std::vector<Ubpa::pointf2> curve2;
	std::vector<Ubpa::pointf2> curve3;
	std::vector<Ubpa::pointf2> curve4;
	Ubpa::valf2 scrolling{ 0.f,0.f };
	bool opt_enable_grid{ true };
	bool opt_enable_context_menu{ true };
	bool adding_line{ false };
};

#include "details/CanvasData_AutoRefl.inl"
