#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include <limits> // For numeric_limits
#include <string>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"
#include <gsl/gsl-lite.hpp>

#include "openmc/boundary_condition.h"
#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Surface;

namespace model {
extern std::unordered_map<int, int> surface_map;
extern vector<unique_ptr<Surface>> surfaces;
} // namespace model

//==============================================================================
//! Coordinates for an axis-aligned cuboid that bounds a geometric object.
//==============================================================================

struct BoundingBox {
  double xmin = -INFTY;
  double xmax = INFTY;
  double ymin = -INFTY;
  double ymax = INFTY;
  double zmin = -INFTY;
  double zmax = INFTY;

  inline BoundingBox operator&(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result &= other;
  }

  inline BoundingBox operator|(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  inline BoundingBox& operator&=(const BoundingBox& other)
  {
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  inline BoundingBox& operator|=(const BoundingBox& other)
  {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }
};

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

class Surface {
public:
  int32_t id_;                                          //!< Unique ID
  std::string name_;                                //!< User-defined name
  std::string type_;                                //!< Type of surface 
  std::shared_ptr<BoundaryCondition> bc_ {nullptr}; //!< Boundary condition
  GeometryType geom_type_;   //!< Geometry type indicator (CSG or DAGMC)
  bool surf_source_ {false}; //!< Activate source banking for the surface?

  explicit Surface(pugi::xml_node surf_node);
  Surface();

  virtual ~Surface() {}

  //! Determine which side of a surface a point lies on.
  //! \param r The 3D Cartesian coordinate of a point.
  //! \param u A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! \return true if the point is on the "positive" side of the surface and
  //!   false otherwise.
  bool sense(Position r, Direction u) const;

  //! Determine the direction of a ray reflected from the surface.
  //! \param[in] r The point at which the ray is incident.
  //! \param[in] u Incident direction of the ray
  //! \param[inout] p Pointer to the particle
  //! \return Outgoing direction of the ray
  virtual Direction reflect(Position r, Direction u, Particle* p) const;

  virtual Direction diffuse_reflect(
    Position r, Direction u, uint64_t* seed) const;

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0.  This member
  //! function evaluates that mathematical function.
  //! \param r A 3D Cartesian coordinate.
  virtual double evaluate(Position r) const = 0;

  //! Compute the distance between a point and the surface along a ray.
  //! \param r A 3D Cartesian coordinate.
  //! \param u The direction of the ray.
  //! \param coincident A hint to the code that the given point should lie
  //!   exactly on the surface.
  virtual double distance(Position r, Direction u, bool coincident) const = 0;

  //! Compute the local outward normal direction of the surface.
  //! \param r A 3D Cartesian coordinate.
  //! \return Normal direction
  virtual Direction normal(Position r) const = 0;

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  //! Get the BoundingBox for this surface.
  virtual BoundingBox bounding_box(bool /*pos_side*/) const { return {}; }
  
  //! Get name
  //! \return Surface name
  const std::string& name() const { return name_; }

  //! Set name
  void set_name(const std::string& name) { name_ = name; }
  
  //! Get ID of material
  //! \return ID of material
  int32_t id() const { return id_; }

  //! Assign a unique ID to the material
  //! \param[in] Unique ID to assign. A value of -1 indicates that an ID
  //!   should be automatically assigned.
  void set_id(int32_t id);
  
  //! Get type
  //! \return Surface type
  const std::string& type() const { return type_; }

  //! Set type
  void set_type(const std::string& type) { type_ = type; }
  
  //! Set n-th coefficient
  virtual void set_coeff(int n, double coeff) = 0;
  
  //! Get n-th coefficient
  virtual double get_coeff(int n) const = 0;

protected:
  virtual void to_hdf5_inner(hid_t group_id) const = 0;
  // Private data members
  gsl::index index_;
};

class CSGSurface : public Surface {
public:
  explicit CSGSurface(pugi::xml_node surf_node);
  CSGSurface();
};

//==============================================================================
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane : public CSGSurface {
public:
  SurfaceXPlane() {};
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public CSGSurface {
public:
  SurfaceYPlane() {};
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double y0_;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public CSGSurface {
public:
  SurfaceZPlane() {};
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double z0_;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public CSGSurface {
public:
  SurfacePlane() {};
  explicit SurfacePlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double A_, B_, C_, D_;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public CSGSurface {
public:
  SurfaceXCylinder() {};
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double y0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public CSGSurface {
public:
  SurfaceYCylinder() {};
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public CSGSurface {
public:
  SurfaceZCylinder() {};
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, radius_;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public CSGSurface {
public:
  SurfaceSphere() {};
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, radius_;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public CSGSurface {
public:
  SurfaceXCone() {};
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public CSGSurface {
public:
  SurfaceYCone() {};
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public CSGSurface {
public:
  SurfaceZCone() {};
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K =
//! 0\f$
//==============================================================================

class SurfaceQuadric : public CSGSurface {
public:
  SurfaceQuadric() {};
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the x direction
//
//! \f$(x-x_0)^2/B^2 + (\sqrt{(y-y_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceXTorus : public CSGSurface {
public:
  SurfaceXTorus() {};
  explicit SurfaceXTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the y direction
//
//! \f$(y-y_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceYTorus : public CSGSurface {
public:
  SurfaceYTorus() {};
  explicit SurfaceYTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the z direction
//
//! \f$(z-z_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (y-y_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceZTorus : public CSGSurface {
public:
  SurfaceZTorus() {};
  explicit SurfaceZTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  void set_coeff(int n, double coeff) override;
  double get_coeff(int n) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
