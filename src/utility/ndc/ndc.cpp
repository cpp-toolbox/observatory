#include "ndc.hpp"

#include <stdexcept>
#include <iostream>
#include <vector>

// Constructor with explicit dimensions and origin
Grid::Grid(int rows, int cols, float width, float height, float origin_x, float origin_y)
    : rows(rows), cols(cols), grid_width(width), grid_height(height), origin_x(origin_x), origin_y(origin_y),
      rect_width(width / cols), rect_height(height / rows) {}

// Constructor using a Rectangle
Grid::Grid(int rows, int cols, const Rectangle &rect)
    : Grid(rows, cols, rect.width, rect.height, rect.center.x, rect.center.y) {}

// Get the rectangle at a specific row and column
Rectangle Grid::get_at(int col, int row) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Index out of range");
    }

    // Calculate the center position of the rectangle relative to the origin
    float x_center = origin_x - grid_width / 2 + rect_width * (col + 0.5f);   // X position
    float y_center = origin_y + grid_height / 2 - rect_height * (row + 0.5f); // Y position

    glm::vec3 center = {x_center, y_center, 0.0f}; // Center in 3D space

    return Rectangle{center, rect_width, rect_height};
}

// Get all rectangles in a specific row
std::vector<Rectangle> Grid::get_row(int row) const {
    if (row < 0 || row >= rows) {
        throw std::out_of_range("Row index out of range");
    }

    std::vector<Rectangle> row_rectangles;
    for (int col = 0; col < cols; ++col) {
        row_rectangles.push_back(get_at(col, row));
    }
    return row_rectangles;
}

// Get all rectangles in a specific column
std::vector<Rectangle> Grid::get_column(int col) const {
    if (col < 0 || col >= cols) {
        throw std::out_of_range("Column index out of range");
    }

    std::vector<Rectangle> col_rectangles;
    for (int row = 0; row < rows; ++row) {
        col_rectangles.push_back(get_at(col, row));
    }
    return col_rectangles;
}
