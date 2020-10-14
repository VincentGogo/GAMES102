#include "CanvasSystem.h"

#include "../Components/CanvasData.h"

#include <_deps/imgui/imgui.h>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;
using namespace Ubpa;

static int CURVE_POINTS = 20000;

void PolyBaseFunction(std::vector<Ubpa::pointf2>& points, std::vector<Ubpa::pointf2>& curve) {
	curve.clear();
	if (points.size() < 2) return;
	int n = points.size();

	// 归一化
	double x_min, x_max, y_min, y_max;
	x_min = x_max = points[0][0];
	y_min = y_max = points[0][1];
	for (int i = 1; i < n; ++i) {
		x_min = min(x_min, (double)points[i][0]);
		x_max = max(x_max, (double)points[i][0]);
		y_min = min(y_min, (double)points[i][1]);
		y_max = max(y_max, (double)points[i][1]);
	}
	double x_range = (x_max - x_min);
	double y_range = (y_max - y_min);

	// Xnn * An1 = Yn1
	Eigen::MatrixXd X(n,n), A(n,1), Y(n,1);

	double x = 0, y = 0;
	for (int i = 0; i < n; ++i) {
		x = (points[i][0] -x_min) / x_range;
		X(i, 0) = 1.0f;
		for (int j = 1; j < n; ++j) {
			X(i, j) = pow(x, j);
		}
		
		Y(i, 0) = (points[i][1] - y_min) / y_range;
	}

	A = X.fullPivLu().solve(Y);

	double step = 1.0 / CURVE_POINTS;
	curve.resize(CURVE_POINTS+1);
	Eigen::MatrixXd X1(1, n), Y1(1, 1);
	for (int i = 0; i <= CURVE_POINTS; ++i) {
		x = i * step;
		X1(0, 0) = 1;
		for (int j = 1; j < n; ++j) {
			X1(0, j) = pow(x, j);
		}
		Y1 = X1 * A;
		y = Y1(0, 0);
		curve[i][0] = x * x_range + x_min;
		curve[i][1] = y * y_range + y_min;
	}

}

void GaussBaseFunction(std::vector<Ubpa::pointf2>& points, std::vector<Ubpa::pointf2>& curve, double sigma = 1.0) {
	curve.clear();
	if (points.size() < 2) return;
	int n = points.size();

	double x_min, x_max, y_min, y_max;
	x_min = x_max = points[0][0];
	y_min = y_max = points[0][1];
	for (int i = 1; i < n; ++i) {
		x_min = min(x_min, (double)points[i][0]);
		x_max = max(x_max, (double)points[i][0]);
		y_min = min(y_min, (double)points[i][1]);
		y_max = max(y_max, (double)points[i][1]);
	}
	double x_range = (x_max - x_min);
	double y_range = (y_max - y_min);

	// G_X(n,n+1) * B(n+1,1) = Y(n,1)
	Eigen::MatrixXd G_X(n, n+1), B(n+1, 1), Y(n, 1);

	double sig2 = 2 * sigma * sigma;
	double x = 0, y = 0;
	for (int i = 0; i < n; ++i) {
		x = (points[i][0] - x_min);
		G_X(i, 0) = 1.0;
		for (int j = 1; j <= n; ++j) {
			double u = (points[j - 1][0] - x_min);
			G_X(i, j) = exp(-(x-u)*(x-u)/sig2);
		}

		Y(i, 0) = (points[i][1] - y_min);
	}

	B = G_X.fullPivLu().solve(Y);

	double step = x_range / CURVE_POINTS;
	curve.resize(CURVE_POINTS+1);
	Eigen::MatrixXd X1(1, n+1), Y1(1, 1);
	for (int i = 0; i <= CURVE_POINTS; ++i) {
		y = B(0, 0); 
		x = i * step;
		X1(0, 0) = 1.0;
		for (int j = 1; j <= n; ++j) {
			double u = (points[j - 1][0] - x_min);
			//y += B(j, 0) * exp(-(x - u) * (x - u) / sig2);
			X1(0, j) = exp(-(x - u) * (x - u) / sig2);
		}
		Y1 = X1 * B;
		y = Y1(0, 0);
		curve[i][0] = x + x_min;
		curve[i][1] = y + y_min;
	}
}

void PolyFittingFunction(std::vector<Ubpa::pointf2>& points, std::vector<Ubpa::pointf2>& curve, int maxPower) {
	curve.clear();
	if (points.size() < 2) return;
	int n = points.size();
	if (maxPower >= n - 1) return PolyBaseFunction(points, curve);
	// 归一化
	double x_min, x_max, y_min, y_max;
	x_min = x_max = points[0][0];
	y_min = y_max = points[0][1];
	for (int i = 1; i < n; ++i) {
		x_min = min(x_min, (double)points[i][0]);
		x_max = max(x_max, (double)points[i][0]);
		y_min = min(y_min, (double)points[i][1]);
		y_max = max(y_max, (double)points[i][1]);
	}
	double x_range = (x_max - x_min);
	double y_range = (y_max - y_min);

	// Xnn * An1 = Yn1
	Eigen::MatrixXd X(n, maxPower+1), A(maxPower+1, 1), Y(n, 1);

	double x = 0, y = 0;
	for (int i = 0; i < n; ++i) {
		x = (points[i][0] - x_min);
		X(i, 0) = 1.0f;
		for (int j = 1; j <= maxPower; ++j) {
			X(i, j) = pow(x, j);
		}

		Y(i, 0) = (points[i][1] - y_min);
	}

	A = (X.transpose() * X).fullPivLu().solve(X.transpose() * Y);
	//A = X.jacobiSvd(ComputeThinU | ComputeThinV).solve(Y);

	double step = x_range / CURVE_POINTS;
	curve.resize(CURVE_POINTS+1);
	Eigen::MatrixXd X1(1, maxPower+1), Y1(1, 1);
	for (int i = 0; i <= CURVE_POINTS; ++i) {
		y = A(0, 0); 
		x = i * step;
		X1(0, 0) = 1.0;
		for (int j = 1; j <= maxPower; ++j) {
			X1(0, j) = pow(x, j);
		}
		Y1 = X1 * A;
		y = Y1(0, 0);
		curve[i][0] = x + x_min;
		curve[i][1] = y + y_min;
	}
}

void RidgeRegression(std::vector<Ubpa::pointf2>& points, std::vector<Ubpa::pointf2>& curve, int maxPower, float lamda) {
	curve.clear();
	if (points.size() < 2) return;
	int n = points.size();
	if (maxPower >= n - 1) return PolyBaseFunction(points, curve);
	// 归一化
	double x_min, x_max, y_min, y_max;
	x_min = x_max = points[0][0];
	y_min = y_max = points[0][1];
	for (int i = 1; i < n; ++i) {
		x_min = min(x_min, (double)points[i][0]);
		x_max = max(x_max, (double)points[i][0]);
		y_min = min(y_min, (double)points[i][1]);
		y_max = max(y_max, (double)points[i][1]);
	}
	double x_range = (x_max - x_min);
	double y_range = (y_max - y_min);

	// Xnn * An1 = Yn1
	Eigen::MatrixXd X(n, maxPower + 1), A(maxPower + 1, 1), Y(n, 1);

	double x = 0, y = 0;
	for (int i = 0; i < n; ++i) {
		x = (points[i][0] - x_min) / x_range;
		X(i, 0) = 1.0f;
		for (int j = 1; j <= maxPower; ++j) {
			X(i, j) = pow(x, j);
		}

		Y(i, 0) = (points[i][1] - y_min) / y_range;
	}

	double l = (double)lamda / x_range;
	auto b = X.transpose() * X + l * MatrixXd::Identity(maxPower + 1, maxPower + 1);
	auto c = X.transpose() * Y;
	A = b.fullPivLu().solve(c);
	//A = b.Inverse() * c;

	double step = 1.0 / CURVE_POINTS;
	curve.resize(CURVE_POINTS+1);
	Eigen::MatrixXd X1(1, maxPower + 1), Y1(1, 1);
	for (int i = 0; i <= CURVE_POINTS; ++i) {
		x = i * step;
		X1(0, 0) = 1.0;
		for (int j = 1; j <= maxPower; ++j) {
			X1(0, j) = pow(x, j);
		}
		Y1 = X1 * A;
		y = Y1(0,0);
		curve[i][0] = x * x_range + x_min;
		curve[i][1] = y * y_range + y_min;
	}
}

void CanvasSystem::OnUpdate(Ubpa::UECS::Schedule& schedule) {
	schedule.RegisterCommand([](Ubpa::UECS::World* w) {
		auto data = w->entityMngr.GetSingleton<CanvasData>();
		if (!data)
			return;

		if (ImGui::Begin("Canvas")) {
			ImGui::Checkbox("Enable grid", &data->opt_enable_grid);
			ImGui::Checkbox("Enable context menu", &data->opt_enable_context_menu);
			static int select = 1;
			ImGui::RadioButton("Polynomial Basis Function (Red Line)", &select, 1);
			ImGui::RadioButton("Gauss Basis Function (Green Line)", &select, 2);
			static float sigma = 1.0f;
			ImGui::SliderFloat("Sigma:", &sigma, 0.1f, 50);
			ImGui::RadioButton("Polynomial Fitting Function (Blue Line)", &select, 3);
			static int maxPower = 4;
			ImGui::SliderInt("MaxPower:", &maxPower, 1, 15);
			ImGui::RadioButton("Ridge Regression (White Line)", &select, 4);
			static float lamda = 0;
			ImGui::SliderFloat("lamda:", &lamda, -50.0f, 50.0f);
			ImGui::Text("Mouse Left: click to add point,\nMouse Right: drag to scroll, click for context menu.");
			ImGui::InputInt("Curve Totoal Points", &CURVE_POINTS, 100, 1000);
			// Typically you would use a BeginChild()/EndChild() pair to benefit from a clipping region + own scrolling.
			// Here we demonstrate that this can be replaced by simple offsetting + custom drawing + PushClipRect/PopClipRect() calls.
			// To use a child window instead we could use, e.g:
			//      ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));      // Disable padding
			//      ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(50, 50, 50, 255));  // Set a background color
			//      ImGui::BeginChild("canvas", ImVec2(0.0f, 0.0f), true, ImGuiWindowFlags_NoMove);
			//      ImGui::PopStyleColor();
			//      ImGui::PopStyleVar();
			//      [...]
			//      ImGui::EndChild();
			if (ImGui::Button("Generate Curve"))
			{
				if (select == 1) {
					PolyBaseFunction(data->points, data->curve1);
				}
				else if (select == 2) {
					GaussBaseFunction(data->points, data->curve2, sigma);
				}
				else if (select == 3) {
					PolyFittingFunction(data->points, data->curve3, maxPower);
				}
				else if (select == 4) {
					RidgeRegression(data->points, data->curve4, maxPower, lamda);
				}
			}
			

			// Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
			ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
			ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
			if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
			if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
			ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

			// Draw border and background color
			ImGuiIO& io = ImGui::GetIO();
			ImDrawList* draw_list = ImGui::GetWindowDrawList();
			draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
			draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

			// This will catch our interactions
			ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
			const bool is_hovered = ImGui::IsItemHovered(); // Hovered
			const bool is_active = ImGui::IsItemActive();   // Held
			const ImVec2 origin(canvas_p0.x + data->scrolling[0], canvas_p0.y + data->scrolling[1]); // Lock scrolled origin
			const pointf2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

			// Add first and second point
			//if (is_hovered && !data->adding_line && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
			//{
			//	data->points.push_back(mouse_pos_in_canvas);
			//	data->points.push_back(mouse_pos_in_canvas);
			//	data->adding_line = true;
			//}
			//if (data->adding_line)
			//{
			//	data->points.back() = mouse_pos_in_canvas;
			//	if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
			//		data->adding_line = false;
			//}

			if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
				data->points.push_back(mouse_pos_in_canvas);
			}

			// Pan (we use a zero mouse threshold when there's no context menu)
			// You may decide to make that threshold dynamic based on whether the mouse is hovering something etc.
			const float mouse_threshold_for_pan = data->opt_enable_context_menu ? -1.0f : 0.0f;
			if (is_active && ImGui::IsMouseDragging(ImGuiMouseButton_Right, mouse_threshold_for_pan))
			{
				data->scrolling[0] += io.MouseDelta.x;
				data->scrolling[1] += io.MouseDelta.y;
			}

			// Context menu (under default mouse threshold)
			ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
			if (data->opt_enable_context_menu && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
				ImGui::OpenPopupContextItem("context");
			if (ImGui::BeginPopup("context"))
			{
				if (data->adding_line)
					data->points.resize(data->points.size() - 2);
				data->adding_line = false;
				if (ImGui::MenuItem("Remove one", NULL, false, data->points.size() > 0)) { data->points.pop_back(); }
				if (ImGui::MenuItem("Remove all", NULL, false, data->points.size() > 0)) {
					data->points.clear(); 
					data->curve1.clear(); 
					data->curve2.clear();
					data->curve3.clear();
					data->curve4.clear();
				}
				ImGui::EndPopup();
			}

			// Draw grid + all lines in the canvas
			draw_list->PushClipRect(canvas_p0, canvas_p1, true);
			if (data->opt_enable_grid)
			{
				const float GRID_STEP = 64.0f;
				for (float x = fmodf(data->scrolling[0], GRID_STEP); x < canvas_sz.x; x += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y), ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
				for (float y = fmodf(data->scrolling[1], GRID_STEP); y < canvas_sz.y; y += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y), ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
			}
			//for (int n = 0; n < data->points.size(); n += 2)
			//	draw_list->AddLine(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n + 1][0], origin.y + data->points[n + 1][1]), IM_COL32(255, 255, 0, 255), 2.0f);
			for (auto &center : data->points) {
				draw_list->AddCircleFilled(ImVec2(origin.x + center[0], origin.y + center[1]), 4.0f, IM_COL32(255, 255, 0, 255));
			}
			
			// draw curve
			for (size_t n = 0; n + 1 < data->curve1.size(); n += 1) {
				draw_list->AddLine(ImVec2(origin.x + data->curve1[n][0], origin.y + data->curve1[n][1]), ImVec2(origin.x + data->curve1[n + 1][0], origin.y + data->curve1[n + 1][1]), IM_COL32(255, 0, 0, 255), 2.0f);
			}

			for (size_t n = 0; n + 1 < data->curve2.size(); n += 1) {
				draw_list->AddLine(ImVec2(origin.x + data->curve2[n][0], origin.y + data->curve2[n][1]), ImVec2(origin.x + data->curve2[n + 1][0], origin.y + data->curve2[n + 1][1]), IM_COL32(0, 255, 0, 255), 2.0f);
			}

			for (size_t n = 0; n + 1 < data->curve3.size(); n += 1) {
				draw_list->AddLine(ImVec2(origin.x + data->curve3[n][0], origin.y + data->curve3[n][1]), ImVec2(origin.x + data->curve3[n + 1][0], origin.y + data->curve3[n + 1][1]), IM_COL32(0, 0, 255, 255), 2.0f);
			}

			for (size_t n = 0; n + 1 < data->curve4.size(); n += 1) {
				draw_list->AddLine(ImVec2(origin.x + data->curve4[n][0], origin.y + data->curve4[n][1]), ImVec2(origin.x + data->curve4[n + 1][0], origin.y + data->curve4[n + 1][1]), IM_COL32(255, 255, 255, 255), 2.0f);
			}
			
			draw_list->PopClipRect();
		}

		ImGui::End();
	});
}
