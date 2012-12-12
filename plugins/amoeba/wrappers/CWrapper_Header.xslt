<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text" omit-xml-declaration="yes"/>

<xsl:variable name="std_namespace_id" select="/GCC_XML/Namespace[@name='std']/@id"/>
<xsl:variable name="openmm_namespace_id" select="/GCC_XML/Namespace[@name='OpenMM']/@id"/>
<xsl:variable name="bool_type_id" select="/GCC_XML/FundamentalType[@name='bool']/@id"/>
<xsl:variable name="double_type_id" select="/GCC_XML/FundamentalType[@name='double']/@id"/>
<xsl:variable name="string_type_id" select="/GCC_XML/*[@name='string' and @context=$std_namespace_id]/@id"/>
<xsl:variable name="const_string_type_id" select="/GCC_XML/CvQualifiedType[@type=$string_type_id]/@id"/>
<xsl:variable name="const_ref_string_type_id" select="/GCC_XML/ReferenceType[@type=$const_string_type_id]/@id"/>
<xsl:variable name="vector_string_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::basic_string')]/@id"/>
<xsl:variable name="vector_vec3_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;OpenMM::Vec3')]/@id"/>
<xsl:variable name="vector_bond_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::pair&lt;int, int')]/@id"/>
<xsl:variable name="vector_2d_int_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::vector&lt;int')]/@id"/>
<xsl:variable name="map_parameter_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and contains(@name, 'double')]/@id"/>
<xsl:variable name="map_property_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and not(contains(@name, 'double'))]/@id"/>
<xsl:variable name="vector_double_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;double')]/@id"/>
<xsl:variable name="vector_int_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;int')]/@id"/>
<xsl:variable name="vector_tortor_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::vector&lt;std::vector&lt;double')]/@id"/>

<!-- Do not generate functions for the following classes -->
<xsl:variable name="skip_classes" select="('Vec3', 'Context', 'Kernel', 'System', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory', 'ContextImpl', 'OpenMMException', 'Force', 'ForceImpl')"/>
<!-- Do not generate the following functions -->
<xsl:variable name="skip_methods" select="('OpenMM_Context_getState', 'OpenMM_Platform_loadPluginsFromDirectory')"/>
<!-- Suppress any function which references any of the following classes -->
<xsl:variable name="hide_classes" select="('Kernel', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory', 'ContextImpl')"/>

<!-- Main loop over all classes in the OpenMM namespace -->
<xsl:template match="/GCC_XML">
#ifndef AMOEBA_OPENMM_CWRAPPER_H_
#define AMOEBA_OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT_AMOEBA
#define OPENMM_EXPORT_AMOEBA
#endif

/* Global Constants */
 <xsl:for-each select="Variable[@context=$openmm_namespace_id]">
static <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@type"/></xsl:call-template><xsl:value-of select="concat(' OpenMM_', @name, ' = ', number(@init), ';')"/>
 </xsl:for-each>

/* Type Declarations */
 <xsl:for-each select="(Class | Struct)[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
typedef struct OpenMM_<xsl:value-of select="concat(@name, '_struct OpenMM_', @name, ';')"/>
 </xsl:for-each>

typedef struct OpenMM_2D_IntArray_struct OpenMM_2D_IntArray;
typedef struct OpenMM_3D_DoubleArray_struct OpenMM_3D_DoubleArray;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_AMOEBA OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_destroy( OpenMM_3D_DoubleArray* array);

<xsl:call-template name="primitive_array">
 <xsl:with-param name="element_type" select="'int'"/>
 <xsl:with-param name="name" select="'OpenMM_IntArray'"/>
</xsl:call-template>

 <!-- Class members -->
 <xsl:for-each select="Class[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
  <xsl:call-template name="class"/>
 </xsl:for-each>

#if defined(__cplusplus)
}
#endif

#endif /*AMOEBA_OPENMM_CWRAPPER_H_*/
</xsl:template>

<!-- Print out the declarations for a (Primitive)Array type -->
<xsl:template name="primitive_array">
 <xsl:param name="element_type"/>
 <xsl:param name="name"/>
/* <xsl:value-of select="$name"/> */
extern OPENMM_EXPORT_AMOEBA <xsl:value-of select="$name"/>* <xsl:value-of select="$name"/>_create(int size);
extern OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name"/>_destroy(<xsl:value-of select="$name"/>* array);
extern OPENMM_EXPORT_AMOEBA int <xsl:value-of select="$name"/>_getSize(const <xsl:value-of select="$name"/>* array);
extern OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name"/>_resize(<xsl:value-of select="$name"/>* array, int size);
extern OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name"/>_append(<xsl:value-of select="$name"/>* array, <xsl:value-of select="$element_type"/> value);
extern OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name"/>_set(<xsl:value-of select="$name"/>* array, int index, <xsl:value-of select="$element_type"/> value);
extern OPENMM_EXPORT_AMOEBA <xsl:value-of select="concat($element_type, ' ', $name)"/>_get(const <xsl:value-of select="$name"/>* array, int index);
</xsl:template>

<!-- Print out information for a class -->
<xsl:template name="class">
 <xsl:variable name="class_name" select="@name"/>
 <xsl:variable name="class_id" select="@id"/>

/* OpenMM::<xsl:value-of select="concat(@name, '*/')"/>
 <!-- Enumerations -->
 <xsl:for-each select="/GCC_XML/Enumeration[@context=$class_id and @access='public']">
  <xsl:call-template name="enumeration">
   <xsl:with-param name="class_name" select="$class_name"/>
  </xsl:call-template>
 </xsl:for-each>
 <!-- Constructors and destructor -->
 <xsl:if test="not(@abstract=1)">
  <xsl:variable name="constructors" select="/GCC_XML/Constructor[@context=$class_id and @access='public' and not(@artificial='1')]"/>
  <xsl:for-each select="$constructors">
   <xsl:call-template name="constructor">
    <xsl:with-param name="suffix" select="if (position() > 1) then concat('_', position()) else ''"/>
   </xsl:call-template>
  </xsl:for-each>
 </xsl:if>
extern OPENMM_EXPORT_AMOEBA void OpenMM_<xsl:value-of select="concat(@name, '_destroy(OpenMM_', @name, '* target);')"/>
 <!-- Methods -->
 <xsl:variable name="methods" select="/GCC_XML/Method[@context=$class_id and @access='public']"/>
 <xsl:for-each select="$methods">
  <xsl:variable name="node" select="."/>
  <!-- The next line is to deal with overloaded methods that have a const and a non-const version. -->
  <xsl:if test="not(@const=1) or empty($methods[@name=$node/@name and not(@const=1)])">
   <xsl:variable name="hide">
    <xsl:call-template name="should_hide"/>
   </xsl:variable>
   <xsl:if test="string-length($hide)=0">
    <xsl:call-template name="method">
     <xsl:with-param name="class_name" select="$class_name"/>
    </xsl:call-template>
   </xsl:if>
  </xsl:if>
 </xsl:for-each>
</xsl:template>

<!-- Print out the declaration for an enumeration -->
<xsl:template name="enumeration">
 <xsl:param name="class_name"/>
typedef enum {
  <xsl:for-each select="EnumValue">
   <xsl:value-of select="concat('OpenMM_', $class_name, '_', @name, ' = ', @init)"/>
   <xsl:if test="position() &lt; last()">, </xsl:if>
  </xsl:for-each>
} OpenMM_<xsl:value-of select="concat($class_name, '_', @name, ';')"/>
</xsl:template>

<!-- Print out the declaration for a constructor -->
<xsl:template name="constructor">
 <xsl:param name="suffix"/>
extern OPENMM_EXPORT_AMOEBA OpenMM_<xsl:value-of select="concat(@name, '* OpenMM_', @name, '_create', $suffix, '(')"/>
  <xsl:for-each select="Argument">
   <xsl:if test="position() > 1">, </xsl:if>
   <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@type"/></xsl:call-template>
   <xsl:value-of select="concat(' ', @name)"/>
  </xsl:for-each>
  <xsl:value-of select="');'"/>
</xsl:template>

<!-- Print out the declaration for a method -->
<xsl:template name="method">
 <xsl:param name="class_name"/>
extern OPENMM_EXPORT_AMOEBA <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@returns"/></xsl:call-template><xsl:value-of select="concat(' OpenMM_', $class_name, '_', @name, '(')"/>
 <xsl:if test="not(@static='1')">
  <xsl:if test="@const='1'">
   <xsl:value-of select="'const '"/>
  </xsl:if>
  <xsl:value-of select="concat('OpenMM_', $class_name, '* target')"/>
 </xsl:if>
 <xsl:variable name="static" select="@static"/>
 <xsl:for-each select="Argument">
  <xsl:if test="position() > 1 or not($static='1')">, </xsl:if>
  <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@type"/></xsl:call-template>
  <xsl:value-of select="concat(' ', @name)"/>
 </xsl:for-each>
 <xsl:value-of select="');'"/>
</xsl:template>

<!-- Print out the description of a type in the wrapper API -->
<xsl:template name="wrap_type">
 <xsl:param name="type_id"/>
 <xsl:variable name="node" select="/GCC_XML/*[@id=$type_id]"/>
 <xsl:choose>
  <xsl:when test="$type_id=$bool_type_id">
   <xsl:value-of select="'OpenMM_Boolean'"/>
  </xsl:when>
  <xsl:when test="$type_id=$string_type_id">
   <xsl:value-of select="'char*'"/>
  </xsl:when>
  <xsl:when test="$type_id=$const_ref_string_type_id">
   <xsl:value-of select="'const char*'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_string_type_id">
   <xsl:value-of select="'OpenMM_StringArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_vec3_type_id">
   <xsl:value-of select="'OpenMM_Vec3Array'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_tortor_type_id">
   <xsl:value-of select="'OpenMM_3D_DoubleArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_bond_type_id">
   <xsl:value-of select="'OpenMM_BondArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_2d_int_type_id">
   <xsl:value-of select="'OpenMM_2D_IntArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$map_parameter_type_id">
   <xsl:value-of select="'OpenMM_ParameterArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$map_property_type_id">
   <xsl:value-of select="'OpenMM_PropertyArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_double_type_id">
   <xsl:value-of select="'OpenMM_DoubleArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_int_type_id">
   <xsl:value-of select="'OpenMM_IntArray'"/>
  </xsl:when>
  <xsl:when test="local-name($node)='ReferenceType' or local-name($node)='PointerType'">
   <xsl:call-template name="wrap_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
   <xsl:value-of select="'*'"/>
  </xsl:when>
  <xsl:when test="local-name($node)='CvQualifiedType'">
   <xsl:if test="$node/@const=1">
    <xsl:value-of select="'const '"/>
   </xsl:if>
   <xsl:call-template name="wrap_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='Enumeration'">
   <xsl:variable name="class_name" select="/GCC_XML/Class[@id=$node/@context]/@name"/>
   <xsl:value-of select="concat('OpenMM_', $class_name, '_', $node/@name)"/>
  </xsl:when>
  <xsl:when test="$node/@context=$openmm_namespace_id">
   <xsl:value-of select="concat('OpenMM_', $node/@name)"/>
  </xsl:when>
  <xsl:otherwise>
   <xsl:value-of select="$node/@name"/>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<!-- Determine whether a method should be hidden -->
<xsl:template name="should_hide">
 <xsl:variable name="class_id" select="@context"/>
 <xsl:variable name="method_name" select="concat('OpenMM_', /GCC_XML/Class[@id=$class_id]/@name, '_', @name)"/>
 <xsl:if test="not(empty(index-of($skip_methods, $method_name)))">1</xsl:if>
 <xsl:call-template name="hide_type">
  <xsl:with-param name="type_id" select="@returns"/>
 </xsl:call-template>
 <xsl:for-each select="Argument">
  <xsl:call-template name="hide_type">
   <xsl:with-param name="type_id" select="@type"/>
  </xsl:call-template>
 </xsl:for-each>
</xsl:template>

<!-- This is called by should_hide.  It generates output if the specified type should be hidden. -->
<xsl:template name="hide_type">
 <xsl:param name="type_id"/>
 <xsl:variable name="node" select="/GCC_XML/*[@id=$type_id]"/>
 <xsl:choose>
  <xsl:when test="local-name($node)='ReferenceType' or local-name($node)='PointerType' or local-name($node)='CvQualifiedType'">
   <xsl:call-template name="hide_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='Enumeration'">
   <xsl:call-template name="hide_type">
    <xsl:with-param name="type_id" select="$node/@context"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="$node/@context=$openmm_namespace_id and not(empty(index-of($hide_classes, $node/@name)))">
   <xsl:value-of select="1"/>
  </xsl:when>
 </xsl:choose>
</xsl:template>

</xsl:stylesheet> 
