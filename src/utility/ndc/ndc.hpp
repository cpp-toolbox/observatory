#ifndef NDC_HPP
#define NDC_HPP

#include <glm/glm.hpp>
#include <vector>
#include "../../graphics/vertex_geometry/vertex_geometry.hpp"

class Grid {
  public:
    /**
     * @brief Constructor to initialize the grid with the given number of rows, columns, width, height, and origin.
     * @param rows Number of rows in the grid.
     * @param cols Number of columns in the grid.
     * @param width Total width of the grid in NDC.
     * @param height Total height of the grid in NDC.
     * @param origin_x X-coordinate of the grid's origin.
     * @param origin_y Y-coordinate of the grid's origin.
     */
    Grid(int rows, int cols, float width = 2.0f, float height = 2.0f, float origin_x = 0.0f, float origin_y = 0.0f);

    /**
     * @brief Constructor to initialize the grid using a Rectangle.
     * @param rows Number of rows in the grid.
     * @param cols Number of columns in the grid.
     * @param rect Rectangle defining the grid's dimensions and position.
     */
    Grid(int rows, int cols, const Rectangle &rect);

    /**
     * @brief Get the rectangle at a specific row and column.
     * @param row Row index (0-based).
     * @param col Column index (0-based).
     * @return The rectangle at the specified row and column.
     * @throws std::out_of_range if row or col is out of bounds.
     */
    Rectangle get_at(int col, int row) const;

    /**
     * @brief Get all rectangles in a specific row.
     * @param row Row index (0-based).
     * @return A vector of rectangles in the specified row.
     * @throws std::out_of_range if the row is out of bounds.
     */
    std::vector<Rectangle> get_row(int row) const;

    /**
     * @brief Get all rectangles in a specific column.
     * @param col Column index (0-based).
     * @return A vector of rectangles in the specified column.
     * @throws std::out_of_range if the column is out of bounds.
     */
    std::vector<Rectangle> get_column(int col) const;

  private:
    int rows;          // Number of rows in the grid
    int cols;          // Number of columns in the grid
    float grid_width;  // Total width of the grid in NDC
    float grid_height; // Total height of the grid in NDC
    float origin_x;    // X-coordinate of the grid's origin
    float origin_y;    // Y-coordinate of the grid's origin
    float rect_width;  // Width of each rectangle
    float rect_height; // Height of each rectangle
};

#endif // NDC_HPP
