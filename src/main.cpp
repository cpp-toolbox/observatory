// aboutodo
// was going to now implement dynamic loading of scripted events so we can test.
// need to get ui working better on the file browser
#include <fmt/core.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/detail/qualifier.hpp>

// TODO document the fact that this has to be in the same file that the implementation header is or doesn't work
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>

#include "graphics/texture_packer_model_loading/texture_packer_model_loading.hpp"
#include "graphics/animated_texture_atlas/animated_texture_atlas.hpp"
#include "graphics/scripted_events/scripted_scene_manager.hpp"
#include "graphics/scripted_transform/scripted_transform.hpp"
#include "graphics/particle_emitter/particle_emitter.hpp"
#include "graphics/vertex_geometry/vertex_geometry.hpp"
#include "graphics/picking_texture/picking_texture.hpp"
#include "graphics/texture_packer/texture_packer.hpp"
#include "graphics/shader_cache/shader_cache.hpp"
#include "graphics/batcher/generated/batcher.hpp"
#include "graphics/fps_camera/fps_camera.hpp"
#include "graphics/cube_map/cube_map.hpp"
#include "graphics/window/window.hpp"
#include "graphics/colors/colors.hpp"
#include "graphics/ui/ui.hpp"

#include "sound_system/sound_system.hpp"

#include "utility/glfw_lambda_callback_manager/glfw_lambda_callback_manager.hpp"
#include "utility/temporal_binary_signal/temporal_binary_signal.hpp"
#include "utility/rigged_model_loading/rigged_model_loading.hpp"
#include "utility/unique_id_generator/unique_id_generator.hpp"
#include "utility/model_loading/model_loading.hpp"
#include "utility/input_state/input_state.hpp"
#include "utility/stopwatch/stopwatch.hpp"
#include "utility/fs_utils/fs_utils.hpp"

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

void draw_ivpntp_object(std::vector<draw_info::IVPNTexturePacked> &packed_object, Transform &object_transform,
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

void draw_ivptp_object(std::vector<draw_info::IVPNTexturePacked> &packed_object, Transform &object_transform,
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

std::vector<int> generate_ui_for_directory(std::filesystem::path &current_directory,
                                           std::filesystem::path &currently_selected_file,
                                           vertex_geometry::Rectangle &main_file_view_rect, UI &filesystem_browser,
                                           TemporalBinarySignal &directory_click_signal,
                                           TemporalBinarySignal &file_click_signal) {

    std::vector<int> doids_for_clickable_textboxes_for_active_directory;

    std::vector<std::filesystem::path> files_and_dirs = list_files_and_directories(current_directory);

    vertex_geometry::Grid file_rows(files_and_dirs.size(), 1, main_file_view_rect);
    auto file_rects = file_rows.get_column(0);

    std::function<void()> on_hover = []() {};

    // note that the hover color doesn't work right now because of the batcher
    for (int i = 0; i < file_rects.size(); i++) {
        std::filesystem::path file_or_directory_path = files_and_dirs.at(i);
        std::function<void()> on_click = [file_or_directory_path, &filesystem_browser, &main_file_view_rect,
                                          &directory_click_signal, &current_directory, &currently_selected_file,
                                          &file_click_signal]() {
            // this if block is important, since we're inside of a lamda which is attached onto a ui element
            // if we tried to erase the ui element we're existing on then we've just killed ourself and the code
            // would not run anymore, thus, we use a signal and then outside this function we react based on the signal
            // in a safer environment
            if (std::filesystem::is_directory(file_or_directory_path)) {
                current_directory = file_or_directory_path;
                directory_click_signal.toggle_state();
            } else { // there are only files and directories, nothing else
                currently_selected_file = file_or_directory_path;
                file_click_signal.toggle_state();
            }
        };
        std::string file_or_directory_path_str = file_or_directory_path.filename().string();
        int oid = filesystem_browser.add_clickable_textbox(on_click, on_hover, file_or_directory_path_str,
                                                           file_rects.at(i), colors.white, colors.grey);
        doids_for_clickable_textboxes_for_active_directory.push_back(oid);
    }
    return doids_for_clickable_textboxes_for_active_directory;
}

struct IVPNTPModel {
    Transform transform;
    std::vector<draw_info::IVPNTexturePacked> packed_model;
};

struct IVPNTPRModelSE {
    // TODO in the future these should probably not be references
    // because the variables would probably die.
    // yeah so this is causing problems now so we need to make this work properly
    RecIvpntRiggedCollector *rirc;
    ScriptedEvent *scripted_event;
    std::vector<IVPNTPRigged> packed_model;
};

/**
 * @class UIRenderSuiteImpl
 * @brief Implementation of the IUIRenderSuite interface.
 */
class UIRenderSuiteImpl : public IUIRenderSuite {
  public:
    Batcher &batcher;

    explicit UIRenderSuiteImpl(Batcher &batcher) : batcher(batcher) {}

    void render_colored_box(const UIRect &cb) override {
        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            cb.id, cb.ivpsc.indices, cb.ivpsc.xyz_positions, cb.ivpsc.rgb_colors,
            cb.modified_signal.has_just_changed());
    }

    void render_text_box(const UITextBox &tb) override {
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
            tb.id, tb.text_drawing_data.indices, tb.text_drawing_data.xyz_positions,
            tb.text_drawing_data.texture_coordinates);

        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            tb.id, tb.background_ivpsc.indices, tb.background_ivpsc.xyz_positions, tb.background_ivpsc.rgb_colors,
            tb.modified_signal.has_just_changed());
    }

    void render_clickable_text_box(const UIClickableTextBox &cr) override {
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
            cr.id, cr.text_drawing_data.indices, cr.text_drawing_data.xyz_positions,
            cr.text_drawing_data.texture_coordinates);

        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            cr.id, cr.ivpsc.indices, cr.ivpsc.xyz_positions, cr.ivpsc.rgb_colors,
            cr.modified_signal.has_just_changed());
    }

    void render_input_box(const UIInputBox &ib) override {
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
            ib.id, ib.text_drawing_data.indices, ib.text_drawing_data.xyz_positions,
            ib.text_drawing_data.texture_coordinates, ib.modified_signal.has_just_changed());

        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            ib.id, ib.background_ivpsc.indices, ib.background_ivpsc.xyz_positions, ib.background_ivpsc.rgb_colors,
            ib.modified_signal.has_just_changed());
    }

    void render_dropdown(const UIDropdown &dd) override {
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
            dd.id, dd.dropdown_text_data.indices, dd.dropdown_text_data.xyz_positions,
            dd.dropdown_text_data.texture_coordinates, dd.modified_signal.has_just_changed());

        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            dd.id, dd.dropdown_background.indices, dd.dropdown_background.xyz_positions,
            dd.dropdown_background.rgb_colors, dd.modified_signal.has_just_changed());
    }

    void render_dropdown_option(const UIDropdown &dd, const draw_info::IVPSolidColor &ivpsc,
                                const draw_info::IVPTextured &ivpt, unsigned int doid) override {
        batcher.transform_v_with_signed_distance_field_text_shader_batcher.queue_draw(
            doid, ivpt.indices, ivpt.xyz_positions, ivpt.texture_coordinates, dd.modified_signal.has_just_changed());

        batcher.absolute_position_with_colored_vertex_shader_batcher.queue_draw(
            doid, ivpsc.indices, ivpsc.xyz_positions, ivpsc.rgb_colors, dd.modified_signal.has_just_changed());
    }
};

void regenerate_texture_packer_and_update_objects(TexturePacker &texture_packer,
                                                  std::vector<std::string> &used_texture_paths,
                                                  ShaderCache &shader_cache,
                                                  TemporalBinarySignal &texture_packer_regen_signal, CubeMap &skybox,
                                                  std::vector<IVPNTPModel> &packed_models) {

    // todo we just need to implement the code which will load from file in repack textures and then
    // update that file and then use that file
    /*add_unique_files(currently_packed_textures_paths, used_texture_paths);*/
    /*repack_textures();*/

    texture_packer.regenerate(used_texture_paths);
    // shader_cache.set_uniform(
    //     ShaderType::
    //         TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
    //     ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES,
    //     texture_packer.texture_index_to_bounding_box);
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
            ivpnt.packed_texture_coordinates =
                texture_packer.get_packed_texture_coordinates(ivpnt.texture_path, ivpnt.original_texture_coordinates);
            ivpnt.packed_texture_index = texture_packer.get_packed_texture_index_of_texture(ivpnt.texture_path);
        }
    }
}

int main() {
    TemporalBinarySignal mouse_clicked_signal;
    InputState input_state;

    // CAMERA STUFF START

    double start_time_ms = 0;
    double end_time_ms = 15 * 1000;
    float catmullrom_tension = 0.5;
    // TODO: have to do late initialization through a pointer so we can delay initializaton, when
    // the scripted transform can be initialized with four points, then this can be corrected
    std::unique_ptr<ScriptedTransform> cmr_scripted_transform = nullptr;
    /*ScriptedTransform cmr_scripted_transform({}, start_time_ms, end_time_ms, catmullrom_tension);*/
    double time_at_which_scripted_camera_was_started = 0;
    bool camera_mode_activated = false;
    bool scripted_camera_running = false;
    std::vector<draw_info::IndexedVertexPositions> camera_keyframes;
    float duration_from_start_to_end_of_spline = 5;
    int selected_camera_keyframe_index = -1;
    // TODO: don't use raw pointer
    draw_info::IndexedVertexPositions *selected_object = nullptr;
    float cam_reach = 3;

    // CAMERA STUFF END

    std::vector<IVPNTPModel> packed_models;
    // as of right now this can only have size 1, because of how the shader works it would require multiple draw calls
    // to render more than one of these so for now just use size 1 to save frames
    std::vector<IVPNTPRModelSE> animated_packed_models_scripted_event;

    bool model_with_animations_open = false;
    // we have this one too because a scripted event can occur without an animation.
    // maybe that makes no sense, need more time to figure that one out either way
    // re-using this also for scripted events
    bool animation_playing = false;

    bool opened_singular_scripted_event = false;

    float timescale = 1;

    // TODO the following "active things" were created as a temporary fix to allow for the creation of
    // objects without having to properly implment rule of 5 on rirc as I'm not yet familiar enough with that yet.
    // remove these later on.
    // additionally I use active_scripted_event for when i want to test a singular scripted event, this is also bad.
    double current_animation_time = 0;
    RecIvpntRiggedCollector active_rirc;
    std::vector<IVPNTRigged> active_ivpntrs;
    std::vector<IVPNTPRigged> active_ivptprs;
    ScriptedEvent active_scripted_event;

    Stopwatch fps_counter;

    unsigned int flame_id = UniqueIDGenerator::generate();
    bool flame_active = false;
    bool cigarette_light_active = false;
    std::vector<glm::vec3> flame_vertices = vertex_geometry::generate_rectangle_vertices(0, 0, 1, 1);
    std::vector<unsigned int> flame_indices = vertex_geometry::generate_rectangle_indices();
    auto flame_normals = vertex_geometry::generate_rectangle_normals();

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::debug);

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("mwe_shader_cache_logs.txt", true);
    file_sink->set_level(spdlog::level::info);

    std::vector<spdlog::sink_ptr> sinks = {console_sink, file_sink};

    Window window;
    bool start_in_fullscreen = true;
    window.initialize_glfw_glad_and_return_window(SCREEN_WIDTH, SCREEN_HEIGHT, "cpp-tbx demo", start_in_fullscreen,
                                                  true, false, true);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    FPSCamera camera(glm::vec3(0, 0, 0), .5, SCREEN_WIDTH, SCREEN_HEIGHT, 90, 0.1, 500);
    std::function<void(unsigned int)> char_callback = [](unsigned int _) {};

    std::function<void(int, int, int, int)> key_callback = [&](int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS || action == GLFW_RELEASE) {
            Key &active_key = *input_state.glfw_code_to_key.at(key);
            bool is_pressed = (action == GLFW_PRESS);
            active_key.pressed_signal.set_signal(is_pressed);
        }

        if (key == GLFW_KEY_M && action == GLFW_PRESS) {
            camera.toggle_mouse_freeze();
            window.toggle_mouse_mode();
        }
    };
    // was about to do moues processing using input state and do that insie of ui rendering stuff
    std::function<void(double, double)> mouse_pos_callback = [&](double xpos, double ypos) {
        camera.mouse_callback(xpos, ypos);
    };
    /*wrap_member_function(camera, &FPSCamera::mouse_callback);*/
    std::function<void(int, int, int)> mouse_button_callback = [&](int button, int action, int mods) {
        std::cout << "mbc" << std::endl;
        if (action == GLFW_PRESS || action == GLFW_RELEASE) {
            std::cout << "mbc inside" << std::endl;
            Key &active_key = *input_state.glfw_code_to_key.at(button);
            std::cout << active_key.string_repr << std::endl;
            bool is_pressed = (action == GLFW_PRESS);
            active_key.pressed_signal.set_signal(is_pressed);
        }

        /*if (button == GLFW_MOUSE_BUTTON_LEFT) {*/
        /*    if (action == GLFW_PRESS) {*/
        /*        mouse_clicked_signal.set_on();*/
        /*    }*/
        /**/
        /*    if (action == GLFW_RELEASE) {*/
        /*        mouse_clicked_signal.set_off();*/
        /*    }*/
        /*}*/
    };

    GLFWLambdaCallbackManager glcm(window.glfw_window, char_callback, key_callback, mouse_pos_callback,
                                   mouse_button_callback);

    std::vector<ShaderType> requested_shaders = {
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
        ShaderType::TRANSFORM_V_WITH_SIGNED_DISTANCE_FIELD_TEXT,
        ShaderType::ABSOLUTE_POSITION_WITH_COLORED_VERTEX,
        ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,
        ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_COLORED_VERTEX,
        ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_OBJECT_ID};

    ShaderCache shader_cache(requested_shaders, sinks);
    Batcher batcher(shader_cache);

    PickingTexture picking_texture;
    picking_texture.initialize(SCREEN_WIDTH, SCREEN_HEIGHT);

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

    const std::filesystem::path textures_directory = "assets";
    std::filesystem::path output_dir = std::filesystem::path("assets") / "packed_textures";
    int container_side_length = 1024;

    std::cout << "before packing" << std::endl;
    TexturePacker texture_packer(textures_directory, output_dir, container_side_length);
    std::cout << "after packing" << std::endl;

    // shader_cache.set_uniform(
    //     ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
    //     ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);

    /*shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,*/
    /*                         ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES,*/
    /*                         texture_packer.texture_index_to_bounding_box);*/

    /*shader_cache.set_uniform(*/
    /*    ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,*/
    /*    ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, texture_packer.texture_index_to_bounding_box);*/

    // here we're setting the texture units to for the packed textures to 1 because thats where the bounding boxes
    // are stored
    //
    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED,
                             ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, 1);

    shader_cache.set_uniform(
        ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
        ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, 1);

    std::cout << "after set bounding boxes" << std::endl;

    TemporalBinarySignal texture_packer_regen_signal;

    // Define the directory path and file extension
    std::filesystem::path cube_map_dir = std::filesystem::path("assets") / "skybox";
    std::string file_extension = "png";

    // Construct the CubeMap using the variables
    CubeMap skybox(cube_map_dir, file_extension, texture_packer);
    std::cout << "after skybox" << std::endl;

    std::string currently_packed_textures_paths = "assets/packed_textures/currently_packed_texture_paths.txt";

    std::filesystem::path font_info_path =
        std::filesystem::path("assets") / "fonts" / "times_64_sdf_atlas_font_info.json";
    std::filesystem::path font_json_path = std::filesystem::path("assets") / "fonts" / "times_64_sdf_atlas.json";
    std::filesystem::path font_image_path = std::filesystem::path("assets") / "fonts" / "times_64_sdf_atlas.png";

    std::cout << "before FONT atlas" << std::endl;

    FontAtlas font_atlas(font_info_path.string(), font_json_path.string(), font_image_path.string(), SCREEN_WIDTH,
                         false, true);

    std::cout << "after FONT atlas" << std::endl;

    /*std::filesystem::path explosion_animation_path =*/
    /*    std::filesystem::path("assets") / "images" / "explosion_animation.png";*/
    /*std::filesystem::path explosion_json_path = std::filesystem::path("assets") / "images" /
     * "explosion_animation.json";*/
    /*AnimatedTextureAtlas explosion_ata(explosion_json_path.string(), explosion_animation_path.string(), 50.0, false,*/
    /*                                   texture_packer);*/

    std::cout << "after explosion ata" << std::endl;

    // START OF UI GENERATION
    UIRenderSuiteImpl ui_render_suite(batcher);

    // rendering functions for UI start

    // rendering functions for UI end

    // TOP BAR START
    vertex_geometry::Grid ui_grid(10, 10);

    UI top_bar(font_atlas);
    auto open_rect = ui_grid.get_at(0, 0);
    auto settings_rect = ui_grid.get_at(1, 0);

    // TODO: this needs to be toggled on when ANIMATION MODE is on
    auto pause_rect = ui_grid.get_at(6, 0);
    auto play_rect = ui_grid.get_at(7, 0);
    auto go_to_start_rect = ui_grid.get_at(8, 0);
    auto time_scale_rect = ui_grid.get_at(9, 0);

    /*auto top_left_rect = create_rectangle_from_top_left(-1, -1, */

    bool file_browser_active = false;

    std::function<void()> on_click = [&]() { file_browser_active = not file_browser_active; };
    std::function<void()> on_hover = []() {};

    top_bar.add_clickable_textbox(on_click, on_hover, "open", open_rect, colors.grey10, colors.darkgreen);

    top_bar.add_textbox("settings", settings_rect, colors.grey10);

    std::function<void()> pause_on_click = [&]() { animation_playing = false; };
    std::function<void()> pause_on_hover = []() {};
    top_bar.add_clickable_textbox(pause_on_click, pause_on_hover, "pause", pause_rect, colors.grey10, colors.darkgreen);

    std::function<void()> play_on_click = [&]() { animation_playing = true; };
    std::function<void()> play_on_hover = []() {};
    top_bar.add_clickable_textbox(play_on_click, play_on_hover, "play", play_rect, colors.grey10, colors.darkgreen);

    std::function<void()> go_to_start_on_click = [&]() {
        current_animation_time = 0;
        for (auto &apmse : animated_packed_models_scripted_event) {
            apmse.scripted_event->reset_processed_state();
        }
        active_scripted_event.reset_processed_state();
        /*explosion_ata.reset_processed_state();*/
    };
    std::function<void()> go_to_start_on_hover = []() {};
    top_bar.add_clickable_textbox(go_to_start_on_click, go_to_start_on_hover, "go to start", go_to_start_rect,
                                  colors.grey10, colors.darkgreen);

    std::function<void(std::string)> time_scale_on_confirm = [&](std::string s) {
        try {
            timescale = std::stof(s);
        } catch (const std::invalid_argument &e) {
            std::cerr << "Invalid input: '" << s << "' is not a valid float.\n";
        } catch (const std::out_of_range &e) {
            std::cerr << "Invalid input: '" << s << "' is out of range for a float.\n";
        }
    };
    top_bar.add_input_box(time_scale_on_confirm, "1", time_scale_rect, colors.grey10, colors.darkgreen);

    // TOP BAR END
    //
    // CAMERA SPLINE MENU START

    UI camera_spline_ui(font_atlas);
    auto top_left = ui_grid.get_at(8, 1);
    auto bottom_right = ui_grid.get_at(9, 2);
    std::vector<vertex_geometry::Rectangle> rects = {top_left, bottom_right};
    auto title_rect = vertex_geometry::get_bounding_rectangle(rects);

    auto tension_label_rect = ui_grid.get_at(8, 3);
    auto tension_input_box_rect = ui_grid.get_at(9, 3);

    auto duration_label_rect = ui_grid.get_at(8, 4);
    auto duration_input_box_rect = ui_grid.get_at(9, 4);

    rects = {ui_grid.get_at(8, 5), ui_grid.get_at(9, 5)};
    auto camera_play_rect = vertex_geometry::get_bounding_rectangle(rects);

    rects = {ui_grid.get_at(8, 6), ui_grid.get_at(9, 6)};
    auto camera_pause_rect = vertex_geometry::get_bounding_rectangle(rects);

    rects = {ui_grid.get_at(8, 7), ui_grid.get_at(9, 7)};
    auto camera_go_to_start_rect = vertex_geometry::get_bounding_rectangle(rects);

    auto playback_percentage = ui_grid.get_at(8, 8);
    auto playback_input_box_rect = ui_grid.get_at(9, 8);

    camera_spline_ui.add_textbox("camera spline", title_rect, colors.grey15);

    std::function<void(std::string)> camera_spline_tension_on_confirm = [&](std::string s) {
        try {
            catmullrom_tension = std::stof(s);
        } catch (const std::invalid_argument &e) {
            std::cerr << "Invalid input: '" << s << "' is not a valid float.\n";
        } catch (const std::out_of_range &e) {
            std::cerr << "Invalid input: '" << s << "' is out of range for a float.\n";
        }
    };

    camera_spline_ui.add_textbox("tension: ", tension_label_rect, colors.grey15);
    camera_spline_ui.add_input_box(camera_spline_tension_on_confirm, "1", tension_input_box_rect, colors.grey15,
                                   colors.darkgreen);

    std::function<void(std::string)> camera_spline_duration_on_confirm = [&](std::string s) {
        try {
            duration_from_start_to_end_of_spline = std::stof(s);
        } catch (const std::invalid_argument &e) {
            std::cerr << "Invalid input: '" << s << "' is not a valid float.\n";
        } catch (const std::out_of_range &e) {
            std::cerr << "Invalid input: '" << s << "' is out of range for a float.\n";
        }
    };

    camera_spline_ui.add_textbox("duration: ", duration_label_rect, colors.grey15);
    camera_spline_ui.add_input_box(camera_spline_duration_on_confirm, "5", duration_input_box_rect, colors.grey15,
                                   colors.darkgreen);

    /*filesystem_browser.add_clickable_textbox(on_click, on_hover, "^", fb.up_a_dir_button, colors.purple,
     * colors.green);*/

    std::function<void()> on_run_click = [&]() {
        scripted_camera_running = true;
        time_at_which_scripted_camera_was_started = glfwGetTime();
    };
    std::function<void()> on_pause_click = [&]() { scripted_camera_running = false; };

    camera_spline_ui.add_clickable_textbox(on_run_click, on_hover, "run", camera_play_rect, colors.grey15,
                                           colors.darkgreen);
    camera_spline_ui.add_clickable_textbox(on_pause_click, on_hover, "pause", camera_pause_rect, colors.grey15,
                                           colors.darkgreen);

    camera_spline_ui.add_textbox("go to start", camera_go_to_start_rect, colors.grey15);

    camera_spline_ui.add_textbox("playback percentage: ", playback_percentage, colors.grey15);
    camera_spline_ui.add_input_box(camera_spline_duration_on_confirm, "1", playback_input_box_rect, colors.grey15,
                                   colors.darkgreen);

    /*camera_spline_ui.add_input_box(camera_spline_tension_on_confirm, "1", tension_rect, colors.grey11,*/
    /*                               colors.darkgreen);*/

    // signals for buttons
    // FILE BROWSER START
    TemporalBinarySignal file_click_signal;
    TemporalBinarySignal directory_click_signal;
    TemporalBinarySignal up_a_dir_signal;

    std::filesystem::path current_directory = get_home_directory();
    std::filesystem::path currently_selected_file = "";

    UI filesystem_browser(font_atlas);
    FileBrowser fb(1.5, 1.5);

    filesystem_browser.add_colored_rectangle(fb.background_rect, colors.gray10);
    int curr_dir_doid =
        filesystem_browser.add_textbox(current_directory.string(), fb.current_directory_rect, colors.gold);
    filesystem_browser.add_colored_rectangle(fb.main_file_view_rect, colors.gray40);
    int selected_file_doid = filesystem_browser.add_textbox("select a file", fb.file_selection_bar, colors.gray40);
    filesystem_browser.add_textbox("x", fb.close_button, colors.darkred);
    /*filesystem_browser.add_textbox("^", up_a_dir_button, colors.purple);*/

    on_click = [&]() {
        up_a_dir_signal.toggle_state();
        current_directory = get_parent_directory(current_directory);
    };
    on_hover = []() {};
    filesystem_browser.add_clickable_textbox(on_click, on_hover, "^", fb.up_a_dir_button, colors.purple, colors.green);

    std::string currrent_directory_str = current_directory.string();
    std::string currently_selected_file_str = currently_selected_file.string();
    std::vector<int> doids_for_clickable_textboxes_for_active_directory =
        generate_ui_for_directory(current_directory, currently_selected_file, fb.main_file_view_rect,
                                  filesystem_browser, directory_click_signal, file_click_signal);

    on_click = [&]() {
        if (currently_selected_file.filename().string() == "scripted_event.json") {
            std::cout << "loaded up scripted event" << std::endl;
            active_scripted_event.load_in_new_scripted_event(currently_selected_file.string());
            opened_singular_scripted_event = true;
        }

        if (has_extension(currently_selected_file, "fbx")) {
            std::filesystem::path target_file_name = "scripted_event.json";

            if (file_exists_in_same_dir(currently_selected_file, target_file_name)) {
                // The target file exists; construct the full path
                std::filesystem::path target_file_path = currently_selected_file.parent_path() / target_file_name;
                std::cout << "File exists: " << target_file_path << std::endl;

                active_ivpntrs = active_rirc.parse_model_into_ivpntrs(currently_selected_file.string());

                std::vector<std::string> used_texture_paths;
                for (auto &ivpnt : active_ivpntrs) {
                    used_texture_paths.push_back(ivpnt.texture_path);
                }

                regenerate_texture_packer_and_update_objects(texture_packer, used_texture_paths, shader_cache,
                                                             texture_packer_regen_signal, skybox, packed_models);

                active_ivptprs = texture_packer_model_loading::convert_ivpnt_to_ivpntpr(active_ivpntrs, texture_packer);
                active_scripted_event.load_in_new_scripted_event(target_file_path.string());

                animated_packed_models_scripted_event.emplace_back(&active_rirc, &active_scripted_event,
                                                                   active_ivptprs);

                model_with_animations_open = true;
            } else {
                std::cout << "File does not exist: " << target_file_name << " in the same directory as "
                          << currently_selected_file << std::endl;
            }
        }
        if (has_extension(currently_selected_file, "obj")) {
            auto model_we_are_loading = model_loading::parse_model_into_ivpnts(currently_selected_file.string(), false);
            Transform crosshair_transform = Transform();

            std::vector<std::string> used_texture_paths;
            for (auto &ivpnt : model_we_are_loading) {
                used_texture_paths.push_back(ivpnt.texture_path);
            }

            regenerate_texture_packer_and_update_objects(texture_packer, used_texture_paths, shader_cache,
                                                         texture_packer_regen_signal, skybox, packed_models);

            // add the new model
            Transform new_model_transform = Transform();
            std::vector<draw_info::IVPNTexturePacked> packed_model =
                texture_packer_model_loading::convert_ivpnt_to_ivpntp(model_we_are_loading, texture_packer);
            IVPNTPModel pm(new_model_transform, packed_model);
            packed_models.push_back(pm);
        }
    };
    on_hover = []() {};

    filesystem_browser.add_clickable_textbox(on_click, on_hover, "open", fb.open_button, colors.darkgreen,
                                             colors.lightgreen);

    // FILE BROWSER END

    // END OF UI GENERATION

    /*AnimatedTextureAtlas animated_texture_atlas("", "assets/images/flame.png", 500.0, texture_packer);*/

    std::cout << "after UI" << std::endl;

    Transform crosshair_transform = Transform();
    crosshair_transform.scale = glm::vec3(.01, .01, .01);

    auto crosshair = model_loading::parse_model_into_ivpnts("assets/crosshair/3d_crosshair.obj", false);
    std::cout << "before crosshair pack" << std::endl;

    std::vector<draw_info::IVPNTexturePacked> packed_crosshair =
        texture_packer_model_loading::convert_ivpnt_to_ivpntp(crosshair, texture_packer);
    std::cout << "after crosshair pack" << std::endl;

    auto lightbulb = model_loading::parse_model_into_ivpnts("assets/lightbulb/lightbulb.obj", false);
    // we have four point lights atm

    std::vector<draw_info::IVPNTexturePacked> packed_lightbulb_1 =
        texture_packer_model_loading::convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<draw_info::IVPNTexturePacked> packed_lightbulb_2 =
        texture_packer_model_loading::convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<draw_info::IVPNTexturePacked> packed_lightbulb_3 =
        texture_packer_model_loading::convert_ivpnt_to_ivpntp(lightbulb, texture_packer);
    std::vector<draw_info::IVPNTexturePacked> packed_lightbulb_4 =
        texture_packer_model_loading::convert_ivpnt_to_ivpntp(lightbulb, texture_packer);

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

    RecIvpntRiggedCollector rirc;
    /*std::vector<IVPNTRigged> smoke_ivpntrs = rirc.parse_model_into_ivpntrs("assets/test/test.fbx");*/
    std::cout << "delme before smokefbx" << std::endl;
    std::vector<IVPNTRigged> smoke_ivpntrs = rirc.parse_model_into_ivpntrs("assets/smoking/smoking.fbx");
    std::cout << "delme after smokefbx" << std::endl;

    std::vector<IVPNTPRigged> smoke_ivptprs =
        texture_packer_model_loading::convert_ivpnt_to_ivpntpr(smoke_ivpntrs, texture_packer);
    std::cout << "delme after smokefbx pack" << std::endl;

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
        // own thing
        {SoundType::EXPLOSION, "assets/sounds/explosion.wav"},

        // reload thing
        {SoundType::BOLT_CLOSE, "assets/sniper_rifle/bolt_close.wav"},
        {SoundType::BOLT_OPEN, "assets/sniper_rifle/bolt_open.wav"},
        {SoundType::BOLT_SLIDE_CLOSE, "assets/sniper_rifle/bolt_slide_close.wav"},
        {SoundType::BOLT_SLIDE_OPEN, "assets/sniper_rifle/bolt_slide_open.wav"},
        {SoundType::IDLE_EXIT, "assets/sniper_rifle/idle_exit.wav"},
        {SoundType::IDLE_RETURN, "assets/sniper_rifle/idle_return.wav"}

    };

    SoundSystem sound_system(100, sound_type_to_file);

    std::cout << "after soundsystem" << std::endl;

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
    ScriptedEvent scripted_event("assets/smoking/scripted_event.json");

    // TODO: you can also used shared pointers here
    IVPNTPRModelSE smoke_animated(&rirc, &scripted_event, smoke_ivptprs);
    /*animated_packed_models_scripted_event.push_back(smoke_animated);*/

    std::cout << "after scripted events" << std::endl;

    std::vector<glm::ivec4> smoke_bone_ids(4, glm::ivec4(0, 0, 0, 0));   // 4 because square
    std::vector<glm::vec4> smoke_bone_weights(4, glm::vec4(0, 0, 0, 0)); // 4 because square

    std::unordered_map<std::string, std::function<void(bool, bool)>> event_callbacks = {
        {"exit_idle",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::IDLE_EXIT, glm::vec3(0.0)); }},
        {"bolt_grab", [&](bool first_call, bool last_call) {}},
        {"bolt_slide_back",
         [&](bool first_call, bool last_call) {
             sound_system.queue_sound(SoundType::BOLT_SLIDE_OPEN, glm::vec3(0.0));
         }},
        {"bolt_rotate_complete",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::BOLT_CLOSE, glm::vec3(0.0)); }},
        {"bolt_rotate",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::BOLT_OPEN, glm::vec3(0.0)); }},
        {"bolt_slide_complete",
         [&](bool first_call, bool last_call) {
             sound_system.queue_sound(SoundType::BOLT_SLIDE_CLOSE, glm::vec3(0.0));
         }},
        {"return_to_idle",
         [&](bool first_call, bool last_call) {
             sound_system.queue_sound(SoundType::IDLE_RETURN, glm::vec3(0.0));
             std::cout << "returning to idle" << std::endl;
         }},

        // explosion stuff
        {"explosion",
         [&](bool first_call, bool last_call) { sound_system.queue_sound(SoundType::EXPLOSION, glm::vec3(0.0)); }},

        // this is all from the smoking stuff
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

    std::cout << "before get packed" << std::endl;

    auto smoke_vertices = vertex_geometry::generate_square_vertices(0, 0, 0.5);
    auto smoke_indices = vertex_geometry::generate_rectangle_indices();
    std::vector<glm::vec2> smoke_local_uvs = vertex_geometry::generate_rectangle_texture_coordinates();

    // Convert paths to native format
    std::string smoke_texture_path = std::filesystem::path("assets/images/smoke_64px.png").make_preferred().string();

    auto smoke_texture_coordinates = texture_packer.get_packed_texture_coordinates(smoke_texture_path, smoke_local_uvs);

    auto smoke_pt_idx = texture_packer.get_packed_texture_index_of_texture(smoke_texture_path);
    std::vector<int> smoke_pt_idxs(4, smoke_pt_idx); // 4 because square

    std::cout << "before while" << std::endl;

    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_TEXTURE_PACKED, ShaderUniformVariable::CAMERA_TO_CLIP,
                             camera.get_projection_matrix());

    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_OBJECT_ID,
                             ShaderUniformVariable::CAMERA_TO_CLIP, camera.get_projection_matrix());
    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_OBJECT_ID,
                             ShaderUniformVariable::WORLD_TO_CAMERA, glm::mat4(1));

    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_COLORED_VERTEX,
                             ShaderUniformVariable::CAMERA_TO_CLIP, camera.get_projection_matrix());
    shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_COLORED_VERTEX,
                             ShaderUniformVariable::WORLD_TO_CAMERA, glm::mat4(1));

    int width, height;

    double previous_time = glfwGetTime();
    while (!glfwWindowShouldClose(window.glfw_window)) {
        fps_counter.press();
        double current_time = glfwGetTime();
        double delta_time = current_time - previous_time;
        previous_time = current_time;

        glfwGetFramebufferSize(window.glfw_window, &width, &height);

        glViewport(0, 0, width, height);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.1, 0.1, 0.1, 1.0);

        // pass uniforms
        if (window.cursor_is_grabbed) {

            camera.process_input(input_state.is_pressed(EKey::LEFT_CONTROL), input_state.is_pressed(EKey::TAB),
                                 input_state.is_pressed(EKey::w), input_state.is_pressed(EKey::a),
                                 input_state.is_pressed(EKey::s), input_state.is_pressed(EKey::d),
                                 input_state.is_pressed(EKey::SPACE), input_state.is_pressed(EKey::LEFT_SHIFT),
                                 delta_time);
        }

        glm::mat4 projection = camera.get_projection_matrix();

        glm::mat4 view;
        glm::mat4 origin_view;
        if (scripted_camera_running and cmr_scripted_transform) {
            view = cmr_scripted_transform->transform.get_transform_matrix();
            double time_since_scripted_camera_started_sec = glfwGetTime() - time_at_which_scripted_camera_was_started;
            double time_since_scripted_camera_started_ms = time_since_scripted_camera_started_sec * 1000;
            cmr_scripted_transform->update(time_since_scripted_camera_started_ms);
            std::cout << cmr_scripted_transform->transform.get_string_repr() << std::endl;

            Transform transform_without_position = cmr_scripted_transform->transform;
            transform_without_position.position = glm::vec3(0);
            origin_view = transform_without_position.get_transform_matrix();
        } else {
            view = camera.get_view_matrix();
            origin_view = camera.get_view_matrix_at(glm::vec3(0));
        }

        /*glm::mat4 origin_view = camera.get_view_matrix_at(glm::vec3(0));*/
        glm::mat4 local_to_world(1.0f);

        // shader_cache.set_uniform(
        //     ShaderType::
        //         TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
        //     ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        // shader_cache.set_uniform(
        //     ShaderType::
        //         TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES_AND_MULTIPLE_LIGHTS,
        //     ShaderUniformVariable::WORLD_TO_CAMERA, view);

        shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_COLORED_VERTEX,
                                 ShaderUniformVariable::WORLD_TO_CAMERA, view);

        shader_cache.set_uniform(ShaderType::CWL_V_TRANSFORMATION_UBOS_1024_WITH_OBJECT_ID,
                                 ShaderUniformVariable::WORLD_TO_CAMERA, view);

        shader_cache.set_uniform(
            ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
            ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        shader_cache.set_uniform(
            ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES,
            ShaderUniformVariable::WORLD_TO_CAMERA, view);

        // sketchy way to get the skybox working
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
        glDepthFunc(GL_LEQUAL); // change depth function so depth test passes when values are equal to depth
                                // buffer's content
        // Array of skybox faces and their corresponding identifiers
        std::vector<std::tuple<int, const draw_info::IVPTexturePacked &>> skybox_faces = {
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
        // set_shader_light_data(camera, shader_cache, false, glm::vec3(0));
        for (auto &pm : packed_models) {
            /*draw_ivpntp_object(pm.packed_model, pm.transform, ltw_matrices, batcher,*/
            /*                   texture_packer_regen_signal.has_just_changed());*/

            draw_ivptp_object(pm.packed_model, pm.transform, ltw_matrices, batcher,
                              texture_packer_regen_signal.has_just_changed());
        }

        if (animation_playing) {
            current_animation_time += timescale * delta_time;
        }

        /*if (opened_singular_scripted_event) {*/
        /*    active_scripted_event.run_scripted_events(current_animation_time, event_callbacks);*/
        /*}*/

        // draw all animated things
        for (auto &apmse : animated_packed_models_scripted_event) {
            // run the scripted event attached to this:

            if (animation_playing) {
                active_scripted_event.run_scripted_events(current_animation_time, event_callbacks);
                /*scripted_event.run_scripted_events(current_animation_time, event_callbacks);*/

                // first we upload the animation matrix
                std::vector<glm::mat4> bone_transformations;
                apmse.rirc->set_bone_transforms(current_animation_time, bone_transformations);

                const unsigned int MAX_BONES_TO_BE_USED = 100;
                ShaderProgramInfo shader_info = shader_cache.get_shader_program(
                    ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES);

                GLint location = glGetUniformLocation(
                    shader_info.id,
                    shader_cache.get_uniform_name(ShaderUniformVariable::BONE_ANIMATION_TRANSFORMS).c_str());

                shader_cache.use_shader_program(
                    ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES);
                glUniformMatrix4fv(location, MAX_BONES_TO_BE_USED, GL_FALSE, glm::value_ptr(bone_transformations[0]));
            }

            // now the model geometry:
            for (auto &ivpntpr : apmse.packed_model) {
                // Populate bone_indices and bone_weights*/
                std::vector<glm::ivec4> bone_indices;
                std::vector<glm::vec4> bone_weights;

                for (const auto &vertex_bone_data : ivpntpr.bone_data) {
                    glm::ivec4 indices(static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[0]),
                                       static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[1]),
                                       static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[2]),
                                       static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[3]));

                    glm::vec4 weights(vertex_bone_data.weight_value_of_this_vertex_wrt_bone[0],
                                      vertex_bone_data.weight_value_of_this_vertex_wrt_bone[1],
                                      vertex_bone_data.weight_value_of_this_vertex_wrt_bone[2],
                                      vertex_bone_data.weight_value_of_this_vertex_wrt_bone[3]);

                    bone_indices.push_back(indices);
                    bone_weights.push_back(weights);
                }

                std::vector<int> packed_texture_indices(ivpntpr.xyz_positions.size(), ivpntpr.packed_texture_index);
                int ptbbi = texture_packer.get_packed_texture_bounding_box_index_of_texture(ivpntpr.texture);
                std::vector<int> packed_texture_bounding_box_indices(ivpntpr.xyz_positions.size(), ptbbi);

                std::vector<unsigned int> ltw_indices(ivpntpr.xyz_positions.size(), ivpntpr.id);
                /*batcher*/
                /*    .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
                /*    .queue_draw(ivpntpr.id, ivpntpr.indices, ltw_indices, bone_indices, bone_weights,
                 * packed_texture_indices, ivpntpr.packed_texture_coordinates, ivpntpr.normals,
                 * ivpntpr.xyz_positions);*/

                batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher
                    .queue_draw(ivpntpr.id, ivpntpr.indices, ltw_indices, bone_indices, bone_weights,
                                packed_texture_indices, ivpntpr.packed_texture_coordinates,
                                packed_texture_bounding_box_indices, ivpntpr.xyz_positions);
            }
        }

        // run scripted events

        /*std::vector<glm::mat4> bone_transformations;*/
        /*float animation_time_sec = glfwGetTime();*/
        /*rirc.set_bone_transforms(animation_time_sec, bone_transformations);*/
        /**/
        /*const unsigned int MAX_BONES_TO_BE_USED = 100;*/
        /*ShaderProgramInfo shader_info = shader_cache.get_shader_program(*/
        /*    ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES);*/
        /**/
        /*GLint location = glGetUniformLocation(*/
        /*    shader_info.id,
         * shader_cache.get_uniform_name(ShaderUniformVariable::BONE_ANIMATION_TRANSFORMS).c_str());*/
        /**/
        /*shader_cache.use_shader_program(*/
        /*    ShaderType::TEXTURE_PACKER_RIGGED_AND_ANIMATED_CWL_V_TRANSFORMATION_UBOS_1024_WITH_TEXTURES);*/
        /*glUniformMatrix4fv(location, MAX_BONES_TO_BE_USED, GL_FALSE, glm::value_ptr(bone_transformations[0]));*/
        /**/
        /*for (auto &ivptr : smoke_ivptprs) {*/
        /*    // Populate bone_indices and bone_weights*/
        /*    std::vector<glm::ivec4> bone_indices;*/
        /*    std::vector<glm::vec4> bone_weights;*/
        /**/
        /*    for (const auto &vertex_bone_data : ivptr.bone_data) {*/
        /*        glm::ivec4
           indices(static_cast<int>(vertex_bone_data.indices_of_bones_that_affect_this_vertex[0]),*/
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
        /*    int ptbbi = texture_packer.get_packed_texture_bounding_box_index_of_texture(ivptr.texture);*/
        /*    std::vector<int> packed_texture_bounding_box_indices(ivptr.xyz_positions.size(), ptbbi);*/
        /**/
        /*    std::vector<unsigned int> ltw_indices(ivptr.xyz_positions.size(), ivptr.id);*/
        /*    /*batcher*/
        /*    /*
           .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/

        /*    /*    .queue_draw(ivptr.id, ivptr.indices, ltw_indices, bone_indices, bone_weights,
           packed_texture_indices,*/
        /*     * ivptr.packed_texture_coordinates, ivptr.normals, ivptr.xyz_positions);*/
        /**/
        /*    batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher*/
        /*        .queue_draw(ivptr.id, ivptr.indices, ltw_indices, bone_indices, bone_weights,
           packed_texture_indices,*/
        /*                    ivptr.packed_texture_coordinates, packed_texture_bounding_box_indices,
           ivptr.xyz_positions);*/
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

        if (animation_playing) {
            /*double ms_curr_time = current_animation_time * 1000;*/
            /**/
            /*std::vector<glm::vec2> packed_tex_coords =*/
            /*    explosion_ata.get_texture_coordinates_of_current_animation_frame(ms_curr_time);*/
            /**/
            /*bool new_coords = false;*/
            /*if (packed_tex_coords != packed_tex_coords_last_tick) {*/
            /*    new_coords = true;*/
            /*    curr_obj_id += 1;*/
            /*}*/
            /**/
            /*packed_tex_coords_last_tick = packed_tex_coords;*/
            /**/
            /*int ptbbi =*/
            /*    texture_packer.get_packed_texture_bounding_box_index_of_texture(explosion_animation_path.string());*/
            /*const std::vector<int> ptbbis(4, ptbbi);*/
            /**/
            /*int pti = texture_packer.get_packed_texture_index_of_texture(explosion_animation_path.string());*/
            /*const std::vector<int> packed_texture_indices(4, pti);*/
            /*std::vector<unsigned int> ltw_mat_idxs(4, 1);*/
            /**/
            /*batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher*/
            /*    .queue_draw(curr_obj_id, flame_indices, ltw_mat_idxs, smoke_bone_ids, smoke_bone_weights,*/
            /*                packed_texture_indices, packed_tex_coords, ptbbis, flame_vertices, true);*/

            /*batcher*/
            /*    .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
            /*    .queue_draw(curr_obj_id, flame_indices, ltw_mat_idxs, smoke_bone_ids, smoke_bone_weights,*/
            /*                packed_texture_indices, packed_tex_coords, ptbbis, flame_normals, flame_vertices, true);*/
        }

        // render ui now

        auto ndc_mouse_pos = camera.mouse.get_ndc_mouse_pos(SCREEN_WIDTH, SCREEN_HEIGHT);
        auto key_strings_just_pressed = input_state.get_just_pressed_key_strings();
        bool delete_action_just_pressed = input_state.is_just_pressed(EKey::DELETE);
        bool confirm_action_just_pressed = input_state.is_just_pressed(EKey::ENTER);
        bool mouse_just_clicked = input_state.is_just_pressed(EKey::LEFT_MOUSE_BUTTON);

        /*void process_and_queue_render_ui(glm::vec2 ndc_mouse_pos, UI &curr_ui, UIRenderSuite &ui_render_suite,*/
        /*                                 const std::vector<std::string> &key_strings_just_pressed,*/
        /*                                 bool delete_action_just_pressed, bool confirm_action_just_pressed,*/
        /*                                 bool mouse_just_clicked, InputState &input_state) {*/

        process_and_queue_render_ui(ndc_mouse_pos, top_bar, ui_render_suite, key_strings_just_pressed,
                                    delete_action_just_pressed, confirm_action_just_pressed, mouse_just_clicked);
        if (file_browser_active) {
            process_and_queue_render_ui(ndc_mouse_pos, filesystem_browser, ui_render_suite,

                                        key_strings_just_pressed, delete_action_just_pressed,
                                        confirm_action_just_pressed, mouse_just_clicked);
        }

        // draw file browser end ^^^^

        // HANDLE fs browser button clicks
        if (directory_click_signal.has_just_changed() or up_a_dir_signal.has_just_changed()) {
            std::string current_directory_str = current_directory.string();
            std::string currently_selected_file_str = currently_selected_file.string();
            filesystem_browser.modify_text_of_a_textbox(curr_dir_doid, current_directory_str);

            // this is a function which updates the thing, that's all
            // remove all the old ui elements
            for (auto doid : doids_for_clickable_textboxes_for_active_directory) {
                filesystem_browser.remove_clickable_textbox(doid);
            }
            doids_for_clickable_textboxes_for_active_directory.clear();
            // create the new ones
            doids_for_clickable_textboxes_for_active_directory =
                generate_ui_for_directory(current_directory, currently_selected_file, fb.main_file_view_rect,
                                          filesystem_browser, directory_click_signal, file_click_signal);
        }

        if (file_click_signal.has_just_changed()) {
            filesystem_browser.modify_text_of_a_textbox(selected_file_doid, currently_selected_file_str);
        }

        // HANDLE fs browser button clicks END

        // CAMERA LOGIC START
        if (input_state.is_just_pressed(EKey::c)) {
            camera_mode_activated = not camera_mode_activated;
        }

        if (camera_mode_activated) {
            process_and_queue_render_ui(ndc_mouse_pos, camera_spline_ui, ui_render_suite,

                                        key_strings_just_pressed, delete_action_just_pressed,
                                        confirm_action_just_pressed, mouse_just_clicked);

            if (input_state.is_just_pressed(EKey::q)) {

                draw_info::IndexedVertexPositions camera_indicator = vertex_geometry::generate_icosphere(1, .1);
                camera_indicator.transform = camera.transform;
                camera_keyframes.push_back(camera_indicator);

                if (cmr_scripted_transform == nullptr) {
                    if (camera_keyframes.size() >= 4) {
                        std::vector<Transform> initial_keyframes;
                        for (const auto &camera_keyframe : camera_keyframes) {
                            initial_keyframes.push_back(camera_keyframe.transform);
                        }
                        cmr_scripted_transform = std::make_unique<ScriptedTransform>(initial_keyframes, start_time_ms,
                                                                                     end_time_ms, catmullrom_tension);
                        std::cout << "init cmr scripted transform" << std::endl;
                    }
                } else {
                    // TODO: add this back
                    /*cmr_scripted_transform->append_keyframe(camera_indicator.transform);*/
                    std::cout << "appended keyframe scripted transform" << std::endl;
                }
            }
        }

        // camera indicators start
        // TODO: eventually don't use this shader make a camera thing in blender and use that as the object.

        int camera_keyframe_indicator_start_index = 500;

        for (int i = 0; i < camera_keyframes.size(); i++) {
            int idx = camera_keyframe_indicator_start_index + i;

            draw_info::IndexedVertexPositions camera_keyframe = camera_keyframes[i];
            ltw_matrices[idx] = camera_keyframe.transform.get_transform_matrix();
            std::vector<unsigned int> camera_keyframe_object_ids(camera_keyframe.xyz_positions.size(), idx);

            std::vector<glm::vec3> cs(camera_keyframe.xyz_positions.size(), colors.black);
            batcher.cwl_v_transformation_ubos_1024_with_colored_vertex_shader_batcher.queue_draw(
                idx, camera_keyframe.indices, camera_keyframe.xyz_positions, cs, camera_keyframe_object_ids);
        }

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        batcher.cwl_v_transformation_ubos_1024_with_colored_vertex_shader_batcher.draw_everything();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // camera indicators end

        // 3d picking start

        picking_texture.enable_writing();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        for (int i = 0; i < camera_keyframes.size(); i++) {
            int idx = camera_keyframe_indicator_start_index + i;

            draw_info::IndexedVertexPositions camera_keyframe = camera_keyframes[i];
            ltw_matrices[idx] = camera_keyframe.transform.get_transform_matrix();
            std::vector<unsigned int> camera_keyframe_object_ids(camera_keyframe.xyz_positions.size(), idx);

            std::vector<glm::vec3> cs(camera_keyframe.xyz_positions.size(), colors.black);
            batcher.cwl_v_transformation_ubos_1024_with_object_id_shader_batcher.queue_draw(
                idx, camera_keyframe.indices, camera_keyframe_object_ids, camera_keyframe.xyz_positions,
                camera_keyframe_object_ids);
        }

        batcher.cwl_v_transformation_ubos_1024_with_object_id_shader_batcher.draw_everything();
        picking_texture.disable_writing();

        if (input_state.is_pressed(EKey::LEFT_MOUSE_BUTTON)) {
            std::cout << "clicked mouse" << std::endl;

            // using center of screen now because we are relying on where you're looking
            unsigned int center_x = SCREEN_WIDTH / 2;
            unsigned int center_y = SCREEN_HEIGHT / 2;
            PickingTexture::PixelInfo clicked_pixel =
                picking_texture.read_pixel(center_x, SCREEN_HEIGHT - 1 - center_y);

            std::cout << clicked_pixel.object_id << std::endl;
            // TODO: the way this index is handled is bad.
            if (clicked_pixel.object_id >= 500 and clicked_pixel.object_id < 1000) {
                selected_camera_keyframe_index = clicked_pixel.object_id - 500;
                // TODO: unsafely using a raw pointer here, may be the source of problems
                // if the item is deleted and we still have a pointer there.
                selected_object = &camera_keyframes[selected_camera_keyframe_index];
            }
        }

        if (input_state.is_pressed(EKey::RIGHT_MOUSE_BUTTON)) {
            selected_object = nullptr;
        }

        if (selected_object != nullptr) {
            selected_object->transform.position =
                camera.transform.position + cam_reach * camera.transform.compute_forward_vector();
        }

        // 3d picking end

        // CAMERA LOGIC END

        // render ui now

        batcher.texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_shader_batcher
            .draw_everything();

        /*batcher*/
        /*    .texture_packer_rigged_and_animated_cwl_v_transformation_ubos_1024_with_textures_and_multiple_lights_shader_batcher*/
        /*    .draw_everything();*/

        glDisable(GL_DEPTH_TEST);
        // this is bad, pack it up soon
        glActiveTexture(GL_TEXTURE0);
        font_atlas.texture_atlas.bind_texture();
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
        glfwSwapBuffers(window.glfw_window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window.glfw_window);

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
