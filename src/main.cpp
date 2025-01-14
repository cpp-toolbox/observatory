#include <fmt/core.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/detail/qualifier.hpp>

#include "graphics/animated_texture_atlas/animated_texture_atlas.hpp"
#include "graphics/batcher/generated/batcher.hpp"
#include "graphics/cube_map/cube_map.hpp"
#include "graphics/fps_camera/fps_camera.hpp"
#include "graphics/scripted_events/scripted_scene_manager.hpp"
#include "graphics/vertex_geometry/vertex_geometry.hpp"
#include "graphics/window/window.hpp"
#include "graphics/shader_cache/shader_cache.hpp"
#include "graphics/particle_emitter/particle_emitter.hpp"
#include "graphics/texture_packer/texture_packer.hpp"
#include "graphics/texture_packer_model_loading/texture_packer_model_loading.hpp"
#include "graphics/ui/ui.hpp"
#include "graphics/colors/colors.hpp"

#include "sound_system/sound_system.hpp"

#include "utility/glfw_lambda_callback_manager/glfw_lambda_callback_manager.hpp"
#include "utility/model_loading/model_loading.hpp"
#include "utility/rigged_model_loading/rigged_model_loading.hpp"
#include "utility/temporal_binary_signal/temporal_binary_signal.hpp"
#include "utility/unique_id_generator/unique_id_generator.hpp"
#include "utility/stopwatch/stopwatch.hpp"
#include "utility/fs_utils/fs_utils.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

// TODO document the fact that this has to be in the same file that the implementation header is or doesn't work
#include <stb_image.h>
#include <stb_image_write.h>

#include <cstdio>
#include <cstdlib>

#include <functional>
#include <iostream>
#include <random>

#include <set>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

Colors colors;
unsigned int SCREEN_WIDTH = 800;
unsigned int SCREEN_HEIGHT = 800;

#include <atomic>
#include <iostream>

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>

/**
 * @brief Adds new file paths to a file without introducing duplicates.
 *
 * @param file_path The path to the file containing existing file paths.
 * @param new_files A vector of strings representing new file paths to be added.
 */
void add_unique_files(const std::string &file_path, const std::vector<std::string> &new_files) {
    std::unordered_set<std::string> file_set; // To store unique file paths.
    std::ifstream input_file(file_path);      // Open the file for reading.

    // Read existing file paths into the set.
    if (input_file.is_open()) {
        std::string line;
        while (std::getline(input_file, line)) {
            file_set.insert(line);
        }
        input_file.close();
    } else {
        std::cerr << "Could not open file: " << file_path << " for reading.\n";
    }

    // Add new files to the set.
    for (const auto &file : new_files) {
        file_set.insert(file);
    }

    // Write the unique file paths back to the file.
    std::ofstream output_file(file_path, std::ios::trunc); // Open the file for overwriting.
    if (output_file.is_open()) {
        for (const auto &file : file_set) {
            output_file << file << '\n';
        }
        output_file.close();
    } else {
        std::cerr << "Could not open file: " << file_path << " for writing.\n";
    }
}

class BlowingSmokeParticleEmitter {
  public:
    ParticleEmitter particle_emitter;

    BlowingSmokeParticleEmitter(unsigned int max_particles, Transform initial_transform)
        : particle_emitter(life_span_lambda(), initial_velocity_lambda(), velocity_change_lambda(), scaling_lambda(),
                           rotation_lambda(), spawn_delay_lambda(), max_particles, initial_transform) {}

  private:
    static std::function<float()> life_span_lambda() {
        return []() -> float {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> dist(1.0f, 3.0f);
            return dist(rng);
        };
    }

    static std::function<glm::vec3()> initial_velocity_lambda() {
        return []() -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> x_dist(0.05, 0.1f); // minor variance
            std::uniform_real_distribution<float> z_dist(0.3, 0.5f);  // blowing
            std::uniform_real_distribution<float> y_dist(0.5f, 0.1f); // upward push

            // Initial upward push with slight lateral drift
            float dx = x_dist(rng);
            float dy = y_dist(rng); // Main upward velocity
            float dz = z_dist(rng);

            return glm::vec3(dx, dy, dz);
            /*return glm::vec3(0, 0, 0);*/
        };
    }

    static std::function<glm::vec3(float, float)> velocity_change_lambda() {
        return [](float life_percentage, float delta_time) -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.0025f, 0.0025f); // Small lateral variance
            std::uniform_real_distribution<float> vertical_dist(0.05f, 0.1f);

            float accel_x = horizontal_dist(rng);
            float accel_y = vertical_dist(rng);
            accel_y = 0;
            float accel_z = horizontal_dist(rng);

            glm::vec3 smoke_push_down = -glm::vec3(accel_x, accel_y, accel_z) * delta_time;
            return smoke_push_down;
        };
    }

    static std::function<float(float)> scaling_lambda() {
        return [](float life_percentage) -> float { return life_percentage * 0.1; };
    }

    static std::function<float(float)> rotation_lambda() {
        return [](float life_percentage) -> float { return life_percentage / 5.0f; };
    }

    static std::function<float()> spawn_delay_lambda() {
        return []() -> float { return 0.05f; };
    }
};

class CigaretteSmokeParticleEmitter {
  public:
    ParticleEmitter particle_emitter;

    CigaretteSmokeParticleEmitter(unsigned int max_particles, Transform initial_transform)
        : particle_emitter(life_span_lambda(), initial_velocity_lambda(), velocity_change_lambda(), scaling_lambda(),
                           rotation_lambda(), spawn_delay_lambda(), max_particles, initial_transform) {}

  private:
    static std::function<float()> life_span_lambda() {
        return []() -> float {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> dist(1.0f, 3.0f);
            return dist(rng);
        };
    }

    static std::function<glm::vec3()> initial_velocity_lambda() {
        return []() -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.10f, 0.10f); // Small lateral variance
            /*std::uniform_real_distribution<float> upward_dist(0.75f, 0.9f);     // upward push rever back to this*/
            std::uniform_real_distribution<float> upward_dist(0.1f, 0.2f); // upward push

            // Initial upward push with slight lateral drift
            float dx = horizontal_dist(rng);
            float dy = upward_dist(rng); // Main upward velocity
            float dz = horizontal_dist(rng);

            return glm::vec3(dx, dy, dz);
            /*return glm::vec3(0, 0, 0);*/
        };
    }

    static std::function<glm::vec3(float, float)> velocity_change_lambda() {
        return [](float life_percentage, float delta_time) -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.0025f, 0.0025f); // Small lateral variance
            std::uniform_real_distribution<float> vertical_dist(0.05f, 0.1f);

            float accel_x = horizontal_dist(rng);
            float accel_y = vertical_dist(rng);
            accel_y = 0;
            float accel_z = horizontal_dist(rng);

            glm::vec3 smoke_push_down = -glm::vec3(accel_x, accel_y, accel_z) * delta_time;
            return smoke_push_down;
        };
    }

    static std::function<float(float)> scaling_lambda() {
        return [](float life_percentage) -> float { return life_percentage * 0.1; };
    }

    static std::function<float(float)> rotation_lambda() {
        return [](float life_percentage) -> float { return life_percentage / 5.0f; };
    }

    static std::function<float()> spawn_delay_lambda() {
        return []() -> float { return 0.10f; };
    }
};

static void error_callback(int error, const char *description) { fprintf(stderr, "Error: %s\n", description); }

void setVec3(GLint unif_loc, const glm::vec3 &value) { glUniform3fv(unif_loc, 1, &value[0]); }
void setFloat(GLint unif_loc, float value) { glUniform1f(unif_loc, value); }

// Light attribute structure
struct PointLightAttributes {
    glm::vec3 position = glm::vec3(0, 0, 0);
    glm::vec3 ambient = glm::vec3(0, 0, 0);
    glm::vec3 diffuse = glm::vec3(0, 0, 0);
    glm::vec3 specular = glm::vec3(0, 0, 0);
    float constant = 1.0f;
    float linear = 0.09f;
    float quadratic = 0.032f;
};

// NOTE we baked in the specular and diffuse into the lights but in reality this is material based
// need to restructure this later
void set_shader_light_data(FPSCamera &camera, ShaderCache &shader_cache, bool is_flame_active,
                           glm::vec3 flame_light_pos) {
    ShaderProgramInfo shader_info = shader_cache.get_shader_program(
        ShaderType::
            TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS);

    shader_cache.use_shader_program(
        ShaderType::
            TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS);

    GLint location = glGetUniformLocation(shader_info.id, "view_pos");
    if (location != -1) {
        setVec3(location, camera.transform.position);
    } else {
        std::cerr << "Warning: Uniform 'view_pos' not found!" << std::endl;
    }

    // Set directional light
    auto set_dir_light = [&](const glm::vec3 &direction, const glm::vec3 &ambient, const glm::vec3 &diffuse,
                             const glm::vec3 &specular) {
        setVec3(glGetUniformLocation(shader_info.id, "dir_light.direction"), direction);
        setVec3(glGetUniformLocation(shader_info.id, "dir_light.ambient"), ambient);
        setVec3(glGetUniformLocation(shader_info.id, "dir_light.diffuse"), diffuse);
        setVec3(glGetUniformLocation(shader_info.id, "dir_light.specular"), specular);
    };
    set_dir_light({-0.2f, -1.0f, -0.3f}, {0.1f, 0.1f, 0.1f}, {0.8f, 0.8f, 0.8f}, {1.0f, 1.0f, 1.0f});

    // Point light data
    std::vector<PointLightAttributes> point_lights = {
        {flame_light_pos, {0.52f, 0.32f, 0.32f}, {0.1f, 0.1f, 0.1f}, {0.4f, 0.4f, 0.4f}, 8.0f, 8.0f, 8.0f},
        {{.8, -.8, .8}, {0.00f, 0.00f, 0.00f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 1.0f, 0.09f, 0.032f},
        {{.8, .8, -.8}, {0.00f, 0.00f, 0.00f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 1.0f, 0.09f, 0.032f},
        {{-.8, .8, .8}, {0.00f, 0.00f, 0.00f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 1.0f, 0.09f, 0.032f},
    };
    PointLightAttributes empty_light = PointLightAttributes();

    // Enhanced flickering effect
    float current_time = static_cast<float>(glfwGetTime());
    float flicker_factor = (sin(current_time * 7.0f) + sin(current_time * 13.0f)) * 0.5f + 0.5f; // Combined sine waves
    float flame_intensity = is_flame_active ? (0.1f + 0.4f * flicker_factor) : 0.0f;

    for (size_t i = 0; i < point_lights.size(); ++i) {
        PointLightAttributes light_data;
        std::string base = "point_lights[" + std::to_string(i) + "].";

        if (i == 0) {
            if (not is_flame_active) {
                light_data = empty_light;
            } else {
                light_data = point_lights[i];
            }
        } else {
            light_data = point_lights[i];
        }

        setVec3(glGetUniformLocation(shader_info.id, (base + "position").c_str()), light_data.position);
        setVec3(glGetUniformLocation(shader_info.id, (base + "ambient").c_str()), light_data.ambient);
        setVec3(glGetUniformLocation(shader_info.id, (base + "diffuse").c_str()), light_data.diffuse);
        setVec3(glGetUniformLocation(shader_info.id, (base + "specular").c_str()), light_data.specular);
        setFloat(glGetUniformLocation(shader_info.id, (base + "constant").c_str()), light_data.constant);
        setFloat(glGetUniformLocation(shader_info.id, (base + "linear").c_str()), light_data.linear);
        setFloat(glGetUniformLocation(shader_info.id, (base + "quadratic").c_str()), light_data.quadratic);
    }

    // Set spot light
    auto set_spot_light = [&](const glm::vec3 &position, const glm::vec3 &direction, const glm::vec3 &ambient,
                              const glm::vec3 &diffuse, const glm::vec3 &specular, float constant, float linear,
                              float quadratic, float cutoff, float outer_cutoff) {
        setVec3(glGetUniformLocation(shader_info.id, "spot_light.position"), position);
        setVec3(glGetUniformLocation(shader_info.id, "spot_light.direction"), direction);
        setVec3(glGetUniformLocation(shader_info.id, "spot_light.ambient"), ambient);
        setVec3(glGetUniformLocation(shader_info.id, "spot_light.diffuse"), diffuse);
        setVec3(glGetUniformLocation(shader_info.id, "spot_light.specular"), specular);
        setFloat(glGetUniformLocation(shader_info.id, "spot_light.constant"), constant);
        setFloat(glGetUniformLocation(shader_info.id, "spot_light.linear"), linear);
        setFloat(glGetUniformLocation(shader_info.id, "spot_light.quadratic"), quadratic);
        setFloat(glGetUniformLocation(shader_info.id, "spot_light.cut_off"), cutoff);
        setFloat(glGetUniformLocation(shader_info.id, "spot_light.outer_cut_off"), outer_cutoff);
    };

    set_spot_light(camera.transform.position, camera.transform.compute_forward_vector(),
                   /*{0.1f, 0.1f, 0.1f}, // ambient light (soft overall lighting)*/
                   /*{0.5f, 0.5f, 0.5f}, // diffuse light (intense light from the spotlight)*/
                   /*{0.5f, 0.5f, 0.5f}, // specular light (highlighted areas with shininess)*/
                   // temporarily turning this off
                   {0.0f, 0.0f, 0.0f}, // ambient light (soft overall lighting)
                   {0.0f, 0.0f, 0.0f}, // diffuse light (intense light from the spotlight)
                   {0.0f, 0.0f, 0.0f}, // specular light (highlighted areas with shininess)
                   1.0f, 0.09f, 0.032f, glm::cos(glm::radians(12.5f)), glm::cos(glm::radians(15.0f)));
}

// Wrapper that automatically creates a lambda for member functions
template <typename T, typename R, typename... Args> auto wrap_member_function(T &obj, R (T::*f)(Args...)) {
    // Return a std::function that wraps the member function in a lambda
    return std::function<R(Args...)>{[&obj, f](Args &&...args) { return (obj.*f)(std::forward<Args>(args)...); }};
}

void draw_ivpntp_object(std::vector<IVPNTexturePacked> &packed_object, Transform &object_transform,
                        glm::mat4 *ltw_matrices, Batcher &batcher, bool replace) {
    int ltw_mat_idx = packed_object[0].id;
    ltw_matrices[ltw_mat_idx] = object_transform.get_transform_matrix();
    for (auto &ivptp : packed_object) {
        // hopefully the matrix at this index is an identity
        std::vector<unsigned int> ltw_indices(ivptp.xyz_positions.size(), ltw_mat_idx);
        std::vector<int> ptis(ivptp.xyz_positions.size(), ivptp.packed_texture_index);
        std::vector<glm::ivec4> blank_bone_ids(ivptp.xyz_positions.size(), glm::ivec4(0, 0, 0, 0));
        std::vector<glm::vec4> blank_bone_weights(ivptp.xyz_positions.size(), glm::vec4(0, 0, 0, 0));
        std::vector<int> ptbbi(ivptp.xyz_positions.size(), ivptp.packed_texture_bounding_box_index);

        batcher
            .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher
            .queue_draw(ivptp.id, ivptp.indices, ltw_indices, blank_bone_ids, blank_bone_weights, ptis,
                        ivptp.packed_texture_coordinates, ptbbi, ivptp.normals, ivptp.xyz_positions, replace);
    }
}

void draw_ivptp_object(std::vector<IVPNTexturePacked> &packed_object, Transform &object_transform,
                       glm::mat4 *ltw_matrices, Batcher &batcher, bool replace) {
    int ltw_mat_idx = packed_object[0].id;
    ltw_matrices[ltw_mat_idx] = object_transform.get_transform_matrix();
    for (auto &ivptp : packed_object) {
        // hopefully the matrix at this index is an identity
        std::vector<unsigned int> ltw_indices(ivptp.xyz_positions.size(), ltw_mat_idx);
        std::vector<int> ptis(ivptp.xyz_positions.size(), ivptp.packed_texture_index);
        std::vector<glm::ivec4> blank_bone_ids(ivptp.xyz_positions.size(), glm::ivec4(0, 0, 0, 0));
        std::vector<glm::vec4> blank_bone_weights(ivptp.xyz_positions.size(), glm::vec4(0, 0, 0, 0));
        std::vector<int> ptbbi(ivptp.xyz_positions.size(), ivptp.packed_texture_bounding_box_index);

        batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher
            .queue_draw(ivptp.id, ivptp.indices, ltw_indices, blank_bone_ids, blank_bone_weights, ptis,
                        ivptp.packed_texture_coordinates, ptbbi, ivptp.xyz_positions, replace);
    }
}

// so that you can attach an item to a bone and keep it attached while animations still play
glm::mat4 get_the_transform_to_attach_an_object_to_a_bone(std::string bone_name, Transform &bone_origin_offset,
                                                          RecIvpntRiggedCollector &rirc) {

    int bone_index = rirc.bone_name_to_unique_index[bone_name];
    BoneInfo bone_info = rirc.bone_unique_idx_to_info[bone_index];
    auto the_transform_that_translates_the_origin_to_the_bones_origin =
        glm::inverse(bone_info.local_space_to_bone_space_in_bind_pose_transformation);

    // put it in the right spot, then git it some translation
    the_transform_that_translates_the_origin_to_the_bones_origin =
        bone_origin_offset.get_transform_matrix() * the_transform_that_translates_the_origin_to_the_bones_origin;
    // then animate it which will work because the emitter is relative to the mesh in bind pose now.
    auto animated_transform = bone_info.local_space_animated_transform_upto_this_bone *
                              the_transform_that_translates_the_origin_to_the_bones_origin;

    return animated_transform;
}

void on_directory_clicked() {}

std::vector<int> generate_ui_for_directory(std::string &current_directory, std::string &currently_selected_file,
                                           Rectangle &main_file_view_rect, UI &filesystem_browser,
                                           TemporalBinarySignal &directory_click_signal,
                                           TemporalBinarySignal &file_click_signal) {

    std::vector<int> doids_for_clickable_textboxes_for_active_directory;

    std::vector<std::string> files_and_dirs = list_files_and_directories(current_directory);

    Grid file_rows(files_and_dirs.size(), 1, main_file_view_rect);
    auto file_rects = file_rows.get_column(0);

    std::function<void()> on_hover = []() {};

    // note that the hover color doesn't work right now because of the batcher
    for (int i = 0; i < file_rects.size(); i++) {
        std::string file_or_directory_path = files_and_dirs.at(i);
        std::function<void()> on_click = [file_or_directory_path, &filesystem_browser, &main_file_view_rect,
                                          &directory_click_signal, &current_directory, &currently_selected_file,
                                          &file_click_signal]() {
            // this if block is important, since we're inside of a lamda which is attached onto a ui element
            // if we tried to erase the ui element we're existing on then we've just killed ourself and the code
            // would not run anymore, thus, we use a signal and then outside this function we react based on the signal
            // in a safer environment
            if (is_directory(file_or_directory_path)) {
                current_directory = file_or_directory_path;
                directory_click_signal.toggle_state();
            } else { // there are only files and directories, nothing else
                currently_selected_file = file_or_directory_path;
                file_click_signal.toggle_state();
            }
        };

        int oid = filesystem_browser.add_clickable_textbox(on_click, on_hover, file_or_directory_path, file_rects.at(i),
                                                           colors.white, colors.grey);
        doids_for_clickable_textboxes_for_active_directory.push_back(oid);
    }
    return doids_for_clickable_textboxes_for_active_directory;
}

struct IVPNTPModel {
    Transform transform;
    std::vector<IVPNTexturePacked> packed_model;
};

int main() {
    TemporalBinarySignal mouse_clicked_signal;

    std::vector<IVPNTPModel> packed_models;

    Stopwatch fps_counter;

    unsigned int flame_id = UniqueIDGenerator::generate();
    bool flame_active = false;
    bool cigarette_light_active = false;
    std::vector<glm::vec3> flame_vertices = generate_rectangle_vertices(0, 0, .03, .03);
    std::vector<unsigned int> flame_indices = generate_rectangle_indices();
    auto flame_normals = generate_rectangle_normals();

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::debug);

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("mwe_shader_cache_logs.txt", true);
    file_sink->set_level(spdlog::level::info);

    std::vector<spdlog::sink_ptr> sinks = {console_sink, file_sink};

    LiveInputState live_input_state;

    GLFWwindow *window =
        initialize_glfw_glad_and_return_window(SCREEN_WIDTH, SCREEN_HEIGHT, "cpp-tbx demo", true, true, false, true);

    /*DivploCubeMap skybox("assets/skyboxes/dusk_land", "png", ShaderType::SKYBOX, shader_cache);*/

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    FPSCamera camera(glm::vec3(0, 0, 0), 50, SCREEN_WIDTH, SCREEN_HEIGHT, 90, 0.1, 500);
    std::function<void(unsigned int)> char_callback = [](unsigned int _) {};
    std::function<void(int, int, int, int)> key_callback = [&](int key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_M && action == GLFW_PRESS) {
            toggle_mouse_mode(window);
        }
    };
    std::function<void(double, double)> mouse_pos_callback = wrap_member_function(camera, &FPSCamera::mouse_callback);
    std::function<void(int, int, int)> mouse_button_callback = [&](int button, int action, int mods) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            if (action == GLFW_PRESS) {
                mouse_clicked_signal.set_on();
            }

            if (action == GLFW_RELEASE) {
                mouse_clicked_signal.set_off();
            }
        }
    };

    GLFWLambdaCallbackManager glcm(window, char_callback, key_callback, mouse_pos_callback, mouse_button_callback);

    std::vector<ShaderType> requested_shaders = {
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
        ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT, ShaderType::ABSOLUTE_POSITION_WITH_COLORED_VERTEX,
        ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED};

    ShaderCache shader_cache(requested_shaders, sinks);
    Batcher batcher(shader_cache);

    auto text_color = glm::vec3(0.5, 0.5, 1);
    float char_width = 0.5;
    float edge_transition = 0.1;

    shader_cache.use_shader_program(ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT);

    shader_cache.set_uniform(ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT, ShaderUniformVariable::TRANSFORM,
                             glm::mat4(1.0f));

    shader_cache.set_uniform(ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT, ShaderUniformVariable::RGB_COLOR,
                             text_color);

    shader_cache.set_uniform(ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT,
                             ShaderUniformVariable::CHARACTER_WIDTH, char_width);

    shader_cache.set_uniform(ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT,
                             ShaderUniformVariable::EDGE_TRANSITION_WIDTH, edge_transition);
    shader_cache.stop_using_shader_program();

    const std::filesystem::path textures_directory = "assets/";
    const std::filesystem::path output_dir = "assets/packed_textures";
    int container_side_length = 4096;

    TexturePacker texture_packer(textures_directory, output_dir, container_side_length);

    shader_cache.set_uniform(
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
        ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);

    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,
                             ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES,
                             texture_packer.texture_index_to_bounding_box);

    shader_cache.set_uniform(
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
        ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);

    TemporalBinarySignal texture_packer_regen_signal;

    CubeMap skybox("assets/skybox", "png", texture_packer);

    std::string currently_packed_textures_paths = "assets/packed_textures/currently_packed_texture_paths.txt";

    AnimatedTextureAtlas animated_texture_atlas("", "assets/images/flame.png", 50.0, texture_packer);

    FontAtlas font_atlas("assets/fonts/times_64_sdf_atlas_font_info.json", "assets/fonts/times_64_sdf_atlas.json",
                         "assets/fonts/times_64_sdf_atlas.png", SCREEN_WIDTH, false, true);

    Grid ui_grid(10, 10);
    UI top_bar(font_atlas);
    auto open_rect = ui_grid.get_at(0, 0);
    auto settings_rect = ui_grid.get_at(1, 0);

    /*auto top_left_rect = create_rectangle_from_top_left(-1, -1, */

    bool file_browser_active = false;

    std::function<void()> on_click = [&]() { file_browser_active = not file_browser_active; };
    std::function<void()> on_hover = []() {};

    top_bar.add_clickable_textbox(on_click, on_hover, "open", open_rect, colors.grey10, colors.darkgreen);
    top_bar.add_textbox("settings", settings_rect, colors.grey10);

    // signals for buttons
    TemporalBinarySignal file_click_signal;
    TemporalBinarySignal directory_click_signal;
    TemporalBinarySignal up_a_dir_signal;

    std::string current_directory = get_home_directory();
    std::string currently_selected_file = "";

    UI filesystem_browser(font_atlas);
    float fsb_height = 1.5;
    float fsb_width = 1.5;
    float fsb_to_side_edge_dist = fsb_width / 2.0;
    float fsb_to_top_edge_dist = fsb_height / 2.0;

    Rectangle background_rect = create_rectangle(0, 0, fsb_width, fsb_height);
    Rectangle current_directory_rect = create_rectangle(0, .8 * fsb_to_top_edge_dist, .8 * fsb_width, .1 * fsb_height);
    Rectangle main_file_view_rect = create_rectangle(0, 0 * fsb_height, .7 * fsb_width, .7 * fsb_height);
    Rectangle file_selection_bar = create_rectangle(0, -.4 * fsb_height, .8 * fsb_width, .1 * fsb_height);
    Rectangle open_button = create_rectangle(.4 * fsb_width, -.4 * fsb_height, .1 * fsb_width, .1 * fsb_height);
    Rectangle close_button = create_rectangle(.4 * fsb_width, .4 * fsb_height, .05 * fsb_width, .05 * fsb_height);
    Rectangle up_a_dir_button = create_rectangle(-.4 * fsb_width, .4 * fsb_height, .05 * fsb_width, .05 * fsb_height);

    filesystem_browser.add_colored_rectangle(background_rect, colors.gray10);
    int curr_dir_doid = filesystem_browser.add_textbox(current_directory, current_directory_rect, colors.gold);
    filesystem_browser.add_colored_rectangle(main_file_view_rect, colors.gray40);
    int selected_file_doid = filesystem_browser.add_textbox("select a file", file_selection_bar, colors.gray40);
    filesystem_browser.add_textbox("x", close_button, colors.darkred);
    /*filesystem_browser.add_textbox("^", up_a_dir_button, colors.purple);*/

    on_click = [&]() {
        up_a_dir_signal.toggle_state();
        current_directory = get_parent_directory(current_directory);
    };
    on_hover = []() {};
    filesystem_browser.add_clickable_textbox(on_click, on_hover, "^", up_a_dir_button, colors.purple, colors.green);

    std::vector<int> doids_for_clickable_textboxes_for_active_directory =
        generate_ui_for_directory(current_directory, currently_selected_file, main_file_view_rect, filesystem_browser,
                                  directory_click_signal, file_click_signal);

    on_click = [&]() {
        if (has_extension(currently_selected_file, "obj")) {
            auto model_we_are_loading = parse_model_into_ivpnts(currently_selected_file, false);

            std::vector<std::string> used_texture_paths;
            for (auto &ivpnt : model_we_are_loading) {
                used_texture_paths.push_back(ivpnt.texture_path);
            }
            // todo we just need to implement the code which will load from file in repack textures and then
            // update that file and then use that file
            /*add_unique_files(currently_packed_textures_paths, used_texture_paths);*/
            /*repack_textures();*/

            texture_packer.regenerate(used_texture_paths);
            shader_cache.set_uniform(
                ShaderType::
                    TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
                ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);
            shader_cache.set_uniform(
                ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
                ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);
            shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,
                                     ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES,
                                     texture_packer.texture_index_to_bounding_box);
            texture_packer_regen_signal.toggle_state();
            skybox.regenerate();

            for (auto &pm : packed_models) {
                for (auto &ivpnt : pm.packed_model) {
                    ivpnt.packed_texture_coordinates = texture_packer.get_packed_texture_coordinates(
                        ivpnt.texture_path, ivpnt.original_texture_coordinates);
                    ivpnt.packed_texture_index = texture_packer.get_packed_texture_index_of_texture(ivpnt.texture_path);
                }
            }

            // add the new model
            Transform new_model_transform = Transform();
            std::vector<IVPNTexturePacked> packed_model = convert_ivpnt_to_ivpntp(model_we_are_loading, texture_packer);
            IVPNTPModel pm(new_model_transform, packed_model);
            packed_models.push_back(pm);
        }
    };
    on_hover = []() {};

    filesystem_browser.add_clickable_textbox(on_click, on_hover, "open", open_button, colors.darkgreen,
                                             colors.lightgreen);

    /*AnimatedTextureAtlas animated_texture_atlas("", "assets/images/flame.png", 500.0, texture_packer);*/

    RecIvpntRiggedCollector rirc;

    Transform crosshair_transform = Transform();
    crosshair_transform.scale = glm::vec3(.01, .01, .01);
    auto crosshair = parse_model_into_ivpnts("assets/crosshair/3d_crosshair.obj", false);
    std::vector<IVPNTexturePacked> packed_crosshair = convert_ivpnt_to_ivpntp(crosshair, texture_packer);

    auto lightbulb = parse_model_into_ivpnts("assets/lightbulb/lightbulb.obj", false);
    // we have four point lights atm
    std::vector<IVPNTexturePacked> packed_lightbulb_1 = convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<IVPNTexturePacked> packed_lightbulb_2 = convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<IVPNTexturePacked> packed_lightbulb_3 = convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<IVPNTexturePacked> packed_lightbulb_4 = convert_ivpnt_to_ivpntp(lightbulb, texture_packer);

    /*{{0, 0, 0}, {0.52f, 0.32f, 0.32f}, {0.1f, 0.1f, 0.1f}, {0.4f, 0.4f, 0.4f}, 1.0f, 0.09f, 0.032f},*/
    /*{{2, -2, 2}, {0.02f, 0.02f, 0.02f}, {0.1f, 0.1f, 0.1f}, {0.4f, 0.4f, 0.4f}, 1.0f, 0.09f, 0.032f},*/
    /*{{2, 2, -2}, {0.02f, 0.02f, 0.02f}, {0.1f, 0.1f, 0.1f}, {0.4f, 0.4f, 0.4f}, 1.0f, 0.09f, 0.032f},*/
    /*{{-2, 2, 2}, {0.02f, 0.02f, 0.02f}, {0.1f, 0.1f, 0.1f}, {0.4f, 0.4f, 0.4f}, 1.0f, 0.09f, 0.032f},*/
    auto lightbulb_1_transform = Transform();
    lightbulb_1_transform.position = glm::vec3(0, 0, 0);

    auto lightbulb_2_transform = Transform();
    lightbulb_2_transform.position = glm::vec3(.8, -.8, .8);

    auto lightbulb_3_transform = Transform();
    lightbulb_3_transform.position = glm::vec3(.8, .8, -.8);

    auto lightbulb_4_transform = Transform();
    lightbulb_4_transform.position = glm::vec3(-.8, .8, .8);

    /*std::vector<IVPNTRigged> smoke_ivpntrs = rirc.parse_model_into_ivpntrs("assets/test/test.fbx");*/
    std::vector<IVPNTRigged> smoke_ivpntrs = rirc.parse_model_into_ivpntrs("assets/smoking/smoking.fbx");
    std::vector<IVPNTPRigged> smoke_ivptprs = convert_ivpnt_to_ivpntpr(smoke_ivpntrs, texture_packer);

    glfwSwapInterval(0);

    GLuint ltw_matrices_gl_name;
    glm::mat4 ltw_matrices[1024];

    // initialize all matrices to identity matrices
    for (int i = 0; i < 1024; ++i) {
        ltw_matrices[i] = glm::mat4(1.0f);
    }

    glGenBuffers(1, &ltw_matrices_gl_name);
    glBindBuffer(GL_UNIFORM_BUFFER, ltw_matrices_gl_name);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(ltw_matrices), ltw_matrices, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ltw_matrices_gl_name);

    std::unordered_map<SoundType, std::string> sound_type_to_file = {
        {SoundType::LIGHTER_FAIL, "assets/sounds/lighter_fail.mp3"},
        {SoundType::LIGHTER_SUCCESS, "assets/sounds/lighter_success.mp3"},
        {SoundType::GRAB, "assets/sounds/grab.wav"},
        {SoundType::KNIFE_GRAB, "assets/sounds/knife_grab.wav"},
        {SoundType::STAB, "assets/sounds/stab.wav"},
        {SoundType::CIGARETTE_BURN, "assets/sounds/cigarette_burning.wav"},
        {SoundType::EXHALE, "assets/sounds/exhale.wav"},
        {SoundType::AMBIENT, "assets/sounds/ambient.wav"},
        {SoundType::WOOSH, "assets/sounds/woosh.wav"},
    };

    SoundSystem sound_system(100, sound_type_to_file);

    sound_system.queue_sound(SoundType::AMBIENT, glm::vec3(0, 0, 0));

    std::vector<glm::vec2> packed_tex_coords_last_tick{};
    int curr_obj_id = 1000;
    int flame_obj_id = 1000;

    Transform spe_transform;
    spe_transform.position = glm::vec3(0, 0, 0);
    spe_transform.scale = glm::vec3(.2, .2, .2);
    CigaretteSmokeParticleEmitter cs_pe(300, spe_transform);
    BlowingSmokeParticleEmitter bs_pe(300, spe_transform);
    // turn off at first
    bs_pe.particle_emitter.stop_emitting_particles();
    cs_pe.particle_emitter.stop_emitting_particles();
    ScriptedEvent scripted_event("assets/smoking/smoking_event.json");

    std::vector<glm::ivec4> smoke_bone_ids(4, glm::ivec4(0, 0, 0, 0));   // 4 because square
    std::vector<glm::vec4> smoke_bone_weights(4, glm::vec4(0, 0, 0, 0)); // 4 because square

    std::unordered_map<std::string, std::function<void(bool, bool)>> event_callbacks = {
        {"grab_pack",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::GRAB, glm::vec3(0.0)); }},
        {"grab_lighter",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::GRAB, glm::vec3(0.0)); }},
        {"ligher_flick_fail",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::LIGHTER_FAIL, glm::vec3(0.0)); }},
        {"lighter_flick_success",
         [&](bool first_call, bool last_call) {
             sound_system.queue_sound(SoundType::LIGHTER_SUCCESS, glm::vec3(0.0));
         }},
        {"lighter_flame",
         [&](bool first_call, bool last_call) {
             if (first_call) {
                 flame_active = true;
             }
             if (last_call) {
                 flame_active = false;
             }
         }},
        {"grab_knife",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::KNIFE_GRAB, glm::vec3(0.0)); }},
        {"knife_throw",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::WOOSH, glm::vec3(0.0)); }},
        {"couch_stab",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::STAB, glm::vec3(0.0)); }},
        {"inhale",
         [&](bool first_call, bool last_call) {
             if (first_call) {
                 sound_system.queue_sound(SoundType::CIGARETTE_BURN, glm::vec3(0.0));
                 cigarette_light_active = true;
                 cs_pe.particle_emitter.stop_emitting_particles();
             }

             if (last_call) {
                 cigarette_light_active = false;
                 cs_pe.particle_emitter.resume_emitting_particles();
                 return;
             }
         }},
        {"exhale",
         [&](bool first_call, bool last_call) {
             if (first_call) {
                 sound_system.queue_sound(SoundType::EXHALE, glm::vec3(0.0));
                 bs_pe.particle_emitter.resume_emitting_particles();
             }

             if (last_call) {
                 bs_pe.particle_emitter.stop_emitting_particles();
             }
         }},
    };

    auto smoke_vertices = generate_square_vertices(0, 0, 0.5);
    auto smoke_indices = generate_rectangle_indices();
    std::vector<glm::vec2> smoke_local_uvs = generate_rectangle_texture_coordinates();
    auto smoke_texture_coordinates =
        texture_packer.get_packed_texture_coordinates("assets/images/smoke_64px.png", smoke_local_uvs);
    auto smoke_pt_idx = texture_packer.get_packed_texture_index_of_texture("assets/images/smoke_64px.png");
    std::vector<int> smoke_pt_idxs(4, smoke_pt_idx); // 4 because square

    int width, height;

    double previous_time = glfwGetTime();
    while (!glfwWindowShouldClose(window)) {
        fps_counter.press();
        double current_time = glfwGetTime();
        double delta_time = current_time - previous_time;
        previous_time = current_time;

        glfwGetFramebufferSize(window, &width, &height);

        glViewport(0, 0, width, height);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.1, 0.1, 0.1, 1.0);

        // pass uniforms
        camera.process_input(window, delta_time);

        glm::mat4 projection = camera.get_projection_matrix();
        glm::mat4 view = camera.get_view_matrix();
        glm::mat4 origin_view = camera.get_view_matrix_at(glm::vec3(0));
        glm::mat4 local_to_world(1.0f);

        shader_cache.set_uniform(
            ShaderType::
                TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
            ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        shader_cache.set_uniform(
            ShaderType::
                TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
            ShaderUniformVariable::WORLD_TO_CAMERA, view);

        shader_cache.set_uniform(
            ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
            ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        shader_cache.set_uniform(
            ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
            ShaderUniformVariable::WORLD_TO_CAMERA, view);

        shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED, ShaderUniformVariable::CAMERA_TO_CLIP,
                                 projection);
        shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,
                                 ShaderUniformVariable::WORLD_TO_CAMERA, origin_view);
        shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED, ShaderUniformVariable::LOCAL_TO_WORLD,
                                 local_to_world);

        /*cs_pe.particle_emitter.update(delta_time, projection * view);*/
        /*auto cs_particles = cs_pe.particle_emitter.get_particles_sorted_by_distance();*/
        /**/
        /*bs_pe.particle_emitter.update(delta_time, projection * view);*/
        /*auto bs_particles = bs_pe.particle_emitter.get_particles_sorted_by_distance();*/
        /**/
        /*// VVV CIG*/
        /**/
        /*auto custom_transform = Transform();*/
        /*custom_transform.position = glm::vec3(.05, 0, -.05);*/
        /*auto smoke_emitter_at_cig_tip_transform =*/
        /*    get_the_transform_to_attach_an_object_to_a_bone("cig_root", custom_transform, rirc);*/
        /**/
        /*cs_pe.particle_emitter.transform.set_transform_matrix(smoke_emitter_at_cig_tip_transform);*/
        /*ltw_matrices[0] = smoke_emitter_at_cig_tip_transform * crosshair_transform.get_transform_matrix();*/
        /**/
        /*glm::vec4 cig_light_pos = smoke_emitter_at_cig_tip_transform * glm::vec4(.05, 0, -.05, 1);*/
        /*glm::vec3 cig_light_pos_3d = glm::vec3(cig_light_pos);*/
        /**/
        /*// ^^^ CIG*/
        /**/
        /*// VVV MOUTH*/
        /*custom_transform = Transform();*/
        /*custom_transform.position = glm::vec3(0, 0.02, .08);*/
        /*auto smoke_emitter_at_mouth_transform =*/
        /*    get_the_transform_to_attach_an_object_to_a_bone("head", custom_transform, rirc);*/
        /**/
        /*bs_pe.particle_emitter.transform.set_transform_matrix(smoke_emitter_at_mouth_transform);*/
        /*/*ltw_matrices[1] = smoke_emitter_at_mouth_transform * crosshair_transform.get_transform_matrix();*/
        /*ltw_matrices[1] = smoke_emitter_at_mouth_transform * crosshair_transform.get_transform_matrix();*/
        /*// ^^^ MOUTH*/
        /*//*/
        /*// VVV LIGHTER*/
        /*custom_transform = Transform();*/
        /*custom_transform.position = glm::vec3(-.02, 0, .05);*/
        /*auto lighter_transform =*/
        /*    get_the_transform_to_attach_an_object_to_a_bone("lighter_root", custom_transform, rirc);*/
        /**/
        /*/*ltw_matrices[packed_crosshair[0].id] = lighter_transform * crosshair_transform.get_transform_matrix();*/
        /*/*ltw_matrices[1] = lighter_transform * crosshair_transform.get_transform_matrix();*/
        /*ltw_matrices[1] = lighter_transform;*/
        /**/
        /*glm::vec4 lighter_flame_pos = lighter_transform * glm::vec4(-.02, 0, .02, 1);*/
        /*glm::vec3 lighter_flame_pos_3d = glm::vec3(lighter_flame_pos);*/
        /*// ^^^ LIGHTER*/
        /**/
        /*if (flame_active) {*/
        /*    set_shader_light_data(camera, shader_cache, true, lighter_flame_pos_3d);*/
        /*} else if (cigarette_light_active) {*/
        /*    set_shader_light_data(camera, shader_cache, true, cig_light_pos_3d);*/
        /*} else {*/
        /*}*/

        // draw skybox
        glDepthMask(GL_FALSE);
        glDepthFunc(
            GL_LEQUAL); // change depth function so depth test passes when values are equal to depth buffer's content
        // Array of skybox faces and their corresponding identifiers
        std::vector<std::tuple<int, const IVPTexturePacked &>> skybox_faces = {
            {0, skybox.top_face},  {1, skybox.bottom_face}, {2, skybox.right_face},
            {3, skybox.left_face}, {4, skybox.front_face},  {5, skybox.back_face}};

        for (const auto &[id, face] : skybox_faces) {
            std::vector<int> ptis(face.xyz_positions.size(), face.packed_texture_index);
            std::vector<int> ptbbi(face.xyz_positions.size(), face.packed_texture_bounding_box_index);

            // Call to the actual function
            batcher.cwl_v_transformation_texture_packed_shader_batcher.queue_draw(
                id, face.indices, face.xyz_positions, ptis, face.packed_texture_coordinates, ptbbi,
                texture_packer_regen_signal.has_just_changed());
        }

        batcher.cwl_v_transformation_texture_packed_shader_batcher.draw_everything();
        glDepthFunc(GL_LESS);
        glDepthMask(GL_TRUE);

        // draw all other objects
        set_shader_light_data(camera, shader_cache, false, glm::vec3(0));
        for (auto &pm : packed_models) {
            /*draw_ivpntp_object(pm.packed_model, pm.transform, ltw_matrices, batcher,*/
            /*                   texture_packer_regen_signal.has_just_changed());*/

            draw_ivptp_object(pm.packed_model, pm.transform, ltw_matrices, batcher,
                              texture_packer_regen_signal.has_just_changed());
        }

        // run scripted events

        /*std::vector<glm::mat4> bone_transformations;*/
        /*float animation_time_sec = glfwGetTime();*/
        /*rirc.set_bone_transforms(animation_time_sec, bone_transformations);*/

        /*const unsigned int MAX_BONES_TO_BE_USED = 100;*/
        /*ShaderProgramInfo shader_info = shader_cache.get_shader_program(*/
        /*    ShaderType::*/
        /*        TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS);*/
        /*GLint location = glGetUniformLocation(*/
        /*    shader_info.id,
         * shader_cache.get_uniform_name(ShaderUniformVariable::BONE_ANIMATION_TRANSFORMS).c_str());*/
        /*glUniformMatrix4fv(location, MAX_BONES_TO_BE_USED, GL_FALSE, glm::value_ptr(bone_transformations[0]));*/

        /*for (auto &ivptr : smoke_ivptprs) {*/
        /*    // Populate bone_indices and bone_weights*/
        /*    std::vector<glm::ivec4> bone_indices;*/
        /*    std::vector<glm::vec4> bone_weights;*/
        /**/
        /*    for (const auto &vertex_bone_data : ivptr.bone_data) {*/
        /*        glm::ivec4 indices(static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[0]),*/
        /*                           static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[1]),*/
        /*                           static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[2]),*/
        /*                           static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[3]));*/
        /**/
        /*        glm::vec4 weights(vertex_bone_data.weight_value_of_this_vertex_wrt_bone[0],*/
        /*                          vertex_bone_data.weight_value_of_this_vertex_wrt_bone[1],*/
        /*                          vertex_bone_data.weight_value_of_this_vertex_wrt_bone[2],*/
        /*                          vertex_bone_data.weight_value_of_this_vertex_wrt_bone[3]);*/
        /**/
        /*        bone_indices.push_back(indices);*/
        /*        bone_weights.push_back(weights);*/
        /*    }*/
        /**/
        /*    std::vector<int> packed_texture_indices(ivptr.xyz_positions.size(), ivptr.packed_texture_index);*/
        /**/
        /*    std::vector<unsigned int> ltw_indices(ivptr.xyz_positions.size(), ivptr.id);*/
        /*    batcher*/
        /*        .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*        .queue_draw(ivptr.id, ivptr.indices, ltw_indices, bone_indices, bone_weights,
         * packed_texture_indices,*/
        /*                    ivptr.packed_texture_coordinates, ivptr.normals, ivptr.xyz_positions);*/
        /*}*/

        /*for (size_t i = 0; i < cs_particles.size(); ++i) {*/
        /**/
        /*    auto &curr_particle = cs_particles[i];*/
        /**/
        /*    //  compute the up vector (assuming we want it to be along the y-axis)*/
        /*    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);*/
        /*    glm::vec3 forward = camera.transform.compute_forward_vector();*/
        /**/
        /*    glm::vec3 right = glm::normalize(glm::cross(up, forward));*/
        /**/
        /*    up = glm::normalize(glm::cross(forward, right));*/
        /**/
        /*    // this makes it billboarded*/
        /*    glm::mat4 rotation_matrix = glm::mat4(1.0f);*/
        /*    rotation_matrix[0] = glm::vec4(right, 0.0f);*/
        /*    rotation_matrix[1] = glm::vec4(up, 0.0f);*/
        /*    rotation_matrix[2] = glm::vec4(-forward, 0.0f); // We negate the direction for correct facing*/
        /**/
        /*    // I think this is bad.*/
        /*    glm::mat4 transform = glm::translate(glm::mat4(1.0f), curr_particle.transform.position);*/
        /*    transform *= rotation_matrix;*/
        /*    transform = glm::scale(transform, curr_particle.transform.scale);*/
        /*    transform = glm::scale(transform, curr_particle.emitter_transform.scale);*/
        /**/
        /*    // temporary*/
        /*    ltw_matrices[curr_particle.id] = transform;*/
        /**/
        /*    if (curr_particle.is_alive()) {*/
        /**/
        /*        auto nv = generate_rectangle_vertices_3d(curr_particle.transform.position,*/
        /*                                                 camera.transform.compute_right_vector(),*/
        /*                                                 camera.transform.compute_up_vector(), 1, 1);*/
        /**/
        /*        std::vector<unsigned int> smoke_ltw_mat_idxs(4, curr_particle.id);*/
        /*        batcher*/
        /*            .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*            .queue_draw(curr_particle.id, smoke_indices, smoke_ltw_mat_idxs, smoke_bone_ids,
         * smoke_bone_weights,*/
        /*                        smoke_pt_idxs, smoke_texture_coordinates, flame_normals, smoke_vertices);*/
        /*    }*/
        /*}*/

        /*for (size_t i = 0; i < bs_particles.size(); ++i) {*/
        /**/
        /*    auto &curr_particle = bs_particles[i];*/
        /**/
        /*    //  compute the up vector (assuming we want it to be along the y-axis)*/
        /*    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);*/
        /*    glm::vec3 forward = camera.transform.compute_forward_vector();*/
        /**/
        /*    glm::vec3 right = glm::normalize(glm::cross(up, forward));*/
        /**/
        /*    up = glm::normalize(glm::cross(forward, right));*/
        /**/
        /*    // this makes it billboarded*/
        /*    glm::mat4 rotation_matrix = glm::mat4(1.0f);*/
        /*    rotation_matrix[0] = glm::vec4(right, 0.0f);*/
        /*    rotation_matrix[1] = glm::vec4(up, 0.0f);*/
        /*    rotation_matrix[2] = glm::vec4(-forward, 0.0f); // We negate the direction for correct facing*/
        /**/
        /*    // I think this is bad.*/
        /*    glm::mat4 transform = glm::translate(glm::mat4(1.0f), curr_particle.transform.position);*/
        /*    transform *= rotation_matrix;*/
        /*    transform = glm::scale(transform, curr_particle.transform.scale);*/
        /*    transform = glm::scale(transform, curr_particle.emitter_transform.scale);*/
        /**/
        /*    // temporary*/
        /*    ltw_matrices[curr_particle.id] = transform;*/
        /**/
        /*    if (curr_particle.is_alive()) {*/
        /**/
        /*        auto nv = generate_rectangle_vertices_3d(curr_particle.transform.position,*/
        /*                                                 camera.transform.compute_right_vector(),*/
        /*                                                 camera.transform.compute_up_vector(), 1, 1);*/
        /**/
        /*        std::vector<unsigned int> smoke_ltw_mat_idxs(4, curr_particle.id);*/
        /*        batcher*/
        /*            .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*            .queue_draw(curr_particle.id, smoke_indices, smoke_ltw_mat_idxs, smoke_bone_ids,
         * smoke_bone_weights,*/
        /*                        smoke_pt_idxs, smoke_texture_coordinates, flame_normals, smoke_vertices);*/
        /*    }*/
        /*}*/

        /*if (true) {*/
        /*    double ms_curr_time = glfwGetTime() * 1000;*/
        /**/
        /*    std::vector<glm::vec2> packed_tex_coords =*/
        /**/
        /*        animated_texture_atlas.get_texture_coordinates_of_current_animation_frame(ms_curr_time);*/
        /**/
        /**/
        /*    bool new_coords = false;*/
        /*    if (packed_tex_coords != packed_tex_coords_last_tick) {*/
        /*        new_coords = true;*/
        /*        curr_obj_id += 1;*/
        /*    }*/
        /**/
        /*    packed_tex_coords_last_tick = packed_tex_coords;*/
        /**/
        /*    const std::vector<int> packed_texture_indices(4, 0);*/
        /*    std::vector<unsigned int> ltw_mat_idxs(4, 1);*/
        /**/
        /*    batcher*/
        /*        .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*        .queue_draw(curr_obj_id, flame_indices, ltw_mat_idxs, smoke_bone_ids, smoke_bone_weights,*/
        /*                    packed_texture_indices, packed_tex_coords, flame_normals, flame_vertices);*/
        /*}*/

        // render ui now
        auto ndc_mouse_pos = top_bar.get_ndc_mouse_pos(SCREEN_WIDTH, SCREEN_HEIGHT, camera.mouse.last_mouse_position_x,
                                                       camera.mouse.last_mouse_position_y);

        // draw top bar start VVV
        top_bar.process_mouse_position(ndc_mouse_pos);

        if (mouse_clicked_signal.is_just_on()) {
            top_bar.process_mouse_just_clicked(ndc_mouse_pos);
        }

        for (auto &tb : top_bar.get_text_boxes()) {
            batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
                tb.id, tb.background_ivpsc.indices, tb.background_ivpsc.xyz_positions, tb.background_ivpsc.rgb_colors);

            batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
                tb.id, tb.text_drawing_data.indices, tb.text_drawing_data.xyz_positions,
                tb.text_drawing_data.texture_coordinates);
        }

        for (auto &tb : top_bar.get_clickable_text_boxes()) {
            batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
                tb.id, tb.ivpsc.indices, tb.ivpsc.xyz_positions, tb.ivpsc.rgb_colors,
                tb.modified_signal.has_just_changed());

            if (tb.modified_signal.has_just_changed()) {
                std::cout << "just entered box" << std::endl;
            }

            batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
                tb.id, tb.text_drawing_data.indices, tb.text_drawing_data.xyz_positions,
                tb.text_drawing_data.texture_coordinates);
        }
        // draw top bar end ^^^

        // draw file browser start VVV
        if (file_browser_active) {

            filesystem_browser.process_mouse_position(ndc_mouse_pos);

            if (mouse_clicked_signal.is_just_on()) {
                filesystem_browser.process_mouse_just_clicked(ndc_mouse_pos);
            }

            for (auto &cb : filesystem_browser.get_colored_boxes()) {
                batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
                    cb.id, cb.ivpsc.indices, cb.ivpsc.xyz_positions, cb.ivpsc.rgb_colors);
            }

            for (auto &tb : filesystem_browser.get_text_boxes()) {
                batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
                    tb.id, tb.background_ivpsc.indices, tb.background_ivpsc.xyz_positions,
                    tb.background_ivpsc.rgb_colors);

                batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
                    tb.id, tb.text_drawing_data.indices, tb.text_drawing_data.xyz_positions,
                    tb.text_drawing_data.texture_coordinates, tb.modified_signal.has_just_changed());

                if (tb.modified_signal.has_just_changed()) {
                    std::cout << "changed dir text" << std::endl;
                }
            }

            for (auto &tb : filesystem_browser.get_clickable_text_boxes()) {
                batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
                    tb.id, tb.ivpsc.indices, tb.ivpsc.xyz_positions, tb.ivpsc.rgb_colors,
                    tb.modified_signal.has_just_changed());

                if (tb.modified_signal.has_just_changed()) {
                    std::cout << "just entered box" << std::endl;
                }

                batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
                    tb.id, tb.text_drawing_data.indices, tb.text_drawing_data.xyz_positions,
                    tb.text_drawing_data.texture_coordinates);
            }
        }

        // draw file browser end ^^^^

        if (directory_click_signal.has_just_changed() or up_a_dir_signal.has_just_changed()) {

            filesystem_browser.modify_text_of_a_textbox(curr_dir_doid, current_directory);

            // this is a function which updates the thing, that's all
            // remove all the old ui elements
            for (auto doid : doids_for_clickable_textboxes_for_active_directory) {
                filesystem_browser.remove_clickable_textbox(doid);
            }
            doids_for_clickable_textboxes_for_active_directory.clear();
            // create the new ones
            doids_for_clickable_textboxes_for_active_directory =
                generate_ui_for_directory(current_directory, currently_selected_file, main_file_view_rect,
                                          filesystem_browser, directory_click_signal, file_click_signal);
        }

        if (file_click_signal.has_just_changed()) {
            filesystem_browser.modify_text_of_a_textbox(selected_file_doid, currently_selected_file);
        }

        // render ui now

        /*double curr_time_sec = glfwGetTime();*/
        /*scripted_event.run_scripted_events(curr_time_sec, event_callbacks);*/

        batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher
            .draw_everything();
        /*batcher*/
        /*    .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*    .draw_everything();*/

        glDisable(GL_DEPTH_TEST);
        batcher.absolute_position_with_colored_vertex_shader_batcher.draw_everything();
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.draw_everything();
        glEnable(GL_DEPTH_TEST);

        sound_system.play_all_sounds();

        // load in the matrices
        glBindBuffer(GL_UNIFORM_BUFFER, ltw_matrices_gl_name);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(ltw_matrices), ltw_matrices);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        /*batcher.texture_packer_cwl_v_transformation_ubos_1024_multiple_lights_shader_batcher.draw_everything();*/
        // -------------------

        TemporalBinarySignal::process_all();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
