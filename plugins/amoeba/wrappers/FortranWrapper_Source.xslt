<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text" omit-xml-declaration="yes"/>

<xsl:variable name="std_namespace_id" select="/GCC_XML/Namespace[@name='std']/@id"/>
<xsl:variable name="openmm_namespace_id" select="/GCC_XML/Namespace[@name='OpenMM']/@id"/>
<xsl:variable name="void_type_id" select="/GCC_XML/FundamentalType[@name='void']/@id"/>
<xsl:variable name="bool_type_id" select="/GCC_XML/FundamentalType[@name='bool']/@id"/>
<xsl:variable name="int_type_id" select="/GCC_XML/FundamentalType[@name='int']/@id"/>
<xsl:variable name="double_type_id" select="/GCC_XML/FundamentalType[@name='double']/@id"/>
<xsl:variable name="char_type_id" select="/GCC_XML/FundamentalType[@name='char']/@id"/>
<xsl:variable name="string_type_id" select="/GCC_XML/*[@name='string' and @context=$std_namespace_id]/@id"/>
<xsl:variable name="const_char_type_id" select="/GCC_XML/CvQualifiedType[@type=$char_type_id]/@id"/>
<xsl:variable name="ptr_const_char_type_id" select="/GCC_XML/PointerType[@type=$const_char_type_id]/@id"/>
<xsl:variable name="const_string_type_id" select="/GCC_XML/CvQualifiedType[@type=$string_type_id]/@id"/>
<xsl:variable name="const_ref_string_type_id" select="/GCC_XML/ReferenceType[@type=$const_string_type_id]/@id"/>
<xsl:variable name="vec3_type_id" select="/GCC_XML/Class[@name='Vec3' and @context=$openmm_namespace_id]/@id"/>
<xsl:variable name="vector_string_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::basic_string')]/@id"/>
<xsl:variable name="vector_vec3_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;OpenMM::Vec3')]/@id"/>
<xsl:variable name="vector_bond_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::pair&lt;int, int')]/@id"/>
<xsl:variable name="vector_2d_int_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::vector&lt;int')]/@id"/>
<xsl:variable name="map_parameter_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and contains(@name, 'double')]/@id"/>
<xsl:variable name="map_property_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and not(contains(@name, 'double'))]/@id"/>
<xsl:variable name="vector_double_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;double')]/@id"/>
<xsl:variable name="vector_int_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;int')]/@id"/>
<xsl:variable name="vector_tortor_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::vector&lt;std::vector&lt;double')]/@id"/>
<xsl:variable name="newline">
<xsl:text>
</xsl:text>
</xsl:variable>


<!-- Do not generate functions for the following classes -->
<xsl:variable name="skip_classes" select="('Vec3', 'Context', 'Kernel', 'System', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory', 'ContextImpl', 'OpenMMException', 'Force', 'ForceImpl')"/>
<!-- Do not generate the following functions -->
<xsl:variable name="skip_methods" select="('OpenMM_Context_getState', 'OpenMM_Platform_loadPluginsFromDirectory')"/>
<!-- Suppress any function which references any of the following classes -->
<xsl:variable name="hide_classes" select="('Kernel', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory', 'ContextImpl')"/>

<!-- Main loop over all classes in the OpenMM namespace -->
<xsl:template match="/GCC_XML">
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "../../../wrappers/OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include &lt;cstring&gt;
#include &lt;vector&gt;

using namespace OpenMM;
using namespace std;

/* Utilities for dealing with Fortran's blank-padded strings. */
static void copyAndPadString(char* dest, const char* source, int length) {
    bool reachedEnd = false;
    for (int i = 0; i &lt; length; i++) {
        if (source[i] == 0)
            reachedEnd = true;
        dest[i] = (reachedEnd ? ' ' : source[i]);
    }
}

static string makeString(const char* fsrc, int length) {
    while (length &amp;&amp; fsrc[length-1]==' ')
        --length;
    return string(fsrc, length);
}

extern "C" {

 <!-- Class members -->
 <xsl:for-each select="Class[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
  <xsl:call-template name="class"/>
 </xsl:for-each>

}
</xsl:template>

<!-- Print out the definitions for a (Primitive)Array type -->
<xsl:template name="primitive_array">
 <xsl:param name="element_type"/>
 <xsl:param name="name"/>
 <xsl:variable name="name_lower" select="lower-case($name)"/>
 <xsl:variable name="name_upper" select="upper-case($name)"/>
/* <xsl:value-of select="$name"/> */
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_create_(<xsl:value-of select="$name"/>*&amp; result, const int&amp; size) {
    result = <xsl:value-of select="$name"/>_create(size);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_CREATE(<xsl:value-of select="$name"/>*&amp; result, const int&amp; size) {
    result = <xsl:value-of select="$name"/>_create(size);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_destroy_(<xsl:value-of select="$name"/>*&amp; array) {
    <xsl:value-of select="$name"/>_destroy(array);
    array = 0;
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_DESTROY(<xsl:value-of select="$name"/>*&amp; array) {
    <xsl:value-of select="$name"/>_destroy(array);
    array = 0;
}
OPENMM_EXPORT_AMOEBA int <xsl:value-of select="$name_lower"/>_getsize_(const <xsl:value-of select="$name"/>* const&amp; array) {
    return <xsl:value-of select="$name"/>_getSize(array);
}
OPENMM_EXPORT_AMOEBA int <xsl:value-of select="$name_upper"/>_GETSIZE(const <xsl:value-of select="$name"/>* const&amp; array) {
    return <xsl:value-of select="$name"/>_getSize(array);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_resize_(<xsl:value-of select="$name"/>* const&amp; array, const int&amp; size) {
    <xsl:value-of select="$name"/>_resize(array, size);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_RESIZE(<xsl:value-of select="$name"/>* const&amp; array, const int&amp; size) {
    <xsl:value-of select="$name"/>_resize(array, size);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_append_(<xsl:value-of select="$name"/>* const&amp; array, <xsl:value-of select="$element_type"/> value) {
    <xsl:value-of select="$name"/>_append(array, value);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_APPEND(<xsl:value-of select="$name"/>* const&amp; array, <xsl:value-of select="$element_type"/> value) {
    <xsl:value-of select="$name"/>_append(array, value);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_set_(<xsl:value-of select="$name"/>* const&amp; array, const int&amp; index, <xsl:value-of select="$element_type"/> value) {
    <xsl:value-of select="$name"/>_set(array, index-1, value);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_SET(<xsl:value-of select="$name"/>* const&amp; array, const int&amp; index, <xsl:value-of select="$element_type"/> value) {
    <xsl:value-of select="$name"/>_set(array, index-1, value);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_lower"/>_get_(const <xsl:value-of select="$name"/>* const&amp; array, const int&amp; index, <xsl:value-of select="$element_type"/>&amp; result) {
    result = <xsl:value-of select="$name"/>_get(array, index-1);
}
OPENMM_EXPORT_AMOEBA void <xsl:value-of select="$name_upper"/>_GET(const <xsl:value-of select="$name"/>* const&amp; array, const int&amp; index, <xsl:value-of select="$element_type"/>&amp; result) {
    result = <xsl:value-of select="$name"/>_get(array, index-1);
}
</xsl:template>

<!-- Print out information for a class -->
<xsl:template name="class">
 <xsl:variable name="class_name" select="@name"/>
 <xsl:variable name="class_id" select="@id"/>

/* OpenMM::<xsl:value-of select="concat(@name, '*/', $newline)"/>
 <!-- Constructors and destructor -->
 <xsl:if test="not(@abstract=1)">
  <xsl:variable name="constructors" select="/GCC_XML/Constructor[@context=$class_id and @access='public' and not(@artificial='1')]"/>
  <xsl:for-each select="$constructors">
   <xsl:variable name="suffix" select="if (position() > 1) then concat('_', position()) else ''"/>
   <xsl:variable name="function_name" select="concat('openmm_', $class_name, '_create', $suffix)"/>
   <xsl:call-template name="constructor">
    <xsl:with-param name="function_name" select="lower-case(concat($function_name, '_'))"/>
    <xsl:with-param name="suffix" select="$suffix"/>
   </xsl:call-template>
   <xsl:call-template name="constructor">
    <xsl:with-param name="function_name" select="upper-case($function_name)"/>
    <xsl:with-param name="suffix" select="$suffix"/>
   </xsl:call-template>
  </xsl:for-each>
 </xsl:if>
 <xsl:variable name="function_name" select="concat('openmm_', $class_name, '_destroy')"/>
 <xsl:call-template name="destructor">
  <xsl:with-param name="function_name" select="lower-case(concat($function_name, '_'))"/>
 </xsl:call-template>
 <xsl:call-template name="destructor">
  <xsl:with-param name="function_name" select="upper-case($function_name)"/>
 </xsl:call-template>
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
    <xsl:variable name="function_name" select="concat('openmm_', $class_name, '_', @name)"/>
    <xsl:call-template name="method">
     <xsl:with-param name="class_name" select="$class_name"/>
     <xsl:with-param name="class_id" select="$class_id"/>
     <xsl:with-param name="function_name" select="lower-case(concat($function_name, '_'))"/>
    </xsl:call-template>
    <xsl:call-template name="method">
     <xsl:with-param name="class_name" select="$class_name"/>
     <xsl:with-param name="class_id" select="$class_id"/>
     <xsl:with-param name="function_name" select="upper-case($function_name)"/>
    </xsl:call-template>
   </xsl:if>
  </xsl:if>
 </xsl:for-each>
</xsl:template>

<!-- Print out the definition of a constructor -->
<xsl:template name="constructor">
 <xsl:param name="function_name"/>
 <xsl:param name="suffix"/>
OPENMM_EXPORT_AMOEBA <xsl:value-of select="concat('void ', $function_name, '(OpenMM_', @name, '*&amp; result')"/>
  <!-- Generate the list of arguments -->
  <xsl:for-each select="Argument">
   <xsl:value-of select="', '"/>
   <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@type"/></xsl:call-template>
   <xsl:variable name="is_handle">
    <xsl:call-template name="is_handle_type">
     <xsl:with-param name="type_id" select="@type"/>
    </xsl:call-template>
   </xsl:variable>
   <xsl:if test="string-length($is_handle)>0">&amp;</xsl:if>
   <xsl:variable name="type_id" select="@type"/>
   <xsl:variable name="type_node" select="/GCC_XML/*[@id=$type_id]"/>
   <xsl:if test="not(local-name($type_node)='ReferenceType' or local-name($type_node)='PointerType')"> const&amp;</xsl:if>
   <xsl:value-of select="concat(' ', @name)"/>
  </xsl:for-each>
  <!-- If any argument is a string, include a length argument -->
  <xsl:for-each select="Argument">
   <xsl:if test="@type=$const_ref_string_type_id">
    <xsl:value-of select="concat(', int ', @name, '_length')"/>
   </xsl:if>
  </xsl:for-each>
  <!-- Now the constructor body -->
  <xsl:value-of select="') {'"/>
    result = OpenMM_<xsl:value-of select="concat(@name, '_create', $suffix, '(')"/>
  <xsl:for-each select="Argument">
   <xsl:if test="position() > 1">, </xsl:if>
   <xsl:choose>
    <xsl:when test="@type=$const_ref_string_type_id">
     <xsl:value-of select="concat('makeString(', @name, ', ', @name, '_length).c_str()')"/>
    </xsl:when>
    <xsl:otherwise>
     <xsl:value-of select="@name"/>
    </xsl:otherwise>
   </xsl:choose>
  </xsl:for-each>
  <xsl:value-of select="');'"/>
}
</xsl:template>

<!-- Print out the definition of a destructor -->
<xsl:template name="destructor">
 <xsl:param name="function_name"/>
OPENMM_EXPORT_AMOEBA <xsl:value-of select="concat('void ', $function_name, '(OpenMM_', @name, '*&amp; destroy) {')"/>
    OpenMM_<xsl:value-of select="concat(@name, '_destroy(destroy);')"/>
    destroy = 0;
}
</xsl:template>

<!-- Print out the definition of a method -->
<xsl:template name="method">
 <xsl:param name="class_name"/>
 <xsl:param name="class_id"/>
 <xsl:param name="function_name"/>
 <!-- First the method signature -->
 <xsl:variable name="has_return" select="@returns=$int_type_id or @returns=$double_type_id"/>
 <xsl:variable name="has_return_arg" select="not($has_return or @returns=$void_type_id)"/>
OPENMM_EXPORT_AMOEBA <xsl:if test="$has_return">
  <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@returns"/></xsl:call-template>
 </xsl:if>
 <xsl:if test="not($has_return)">
  <xsl:value-of select="'void'"/>
 </xsl:if>
 <xsl:value-of select="concat(' ', $function_name, '(')"/>
 <xsl:if test="not(@static='1')">
  <xsl:if test="@const='1'">
   <xsl:value-of select="'const '"/>
  </xsl:if>
  <xsl:value-of select="concat('OpenMM_', $class_name, '*&amp; target')"/>
 </xsl:if>
 <!-- Generate the list of arguments -->
 <xsl:variable name="method" select="."/>
 <xsl:for-each select="Argument">
  <xsl:if test="position() > 1 or not($method/@static='1')">, </xsl:if>
  <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@type"/></xsl:call-template>
  <xsl:variable name="is_handle">
   <xsl:call-template name="is_handle_type">
    <xsl:with-param name="type_id" select="@type"/>
   </xsl:call-template>
  </xsl:variable>
  <xsl:if test="string-length($is_handle)>0">&amp;</xsl:if>
  <xsl:variable name="type_id" select="@type"/>
  <xsl:variable name="type_node" select="/GCC_XML/*[@id=$type_id]"/>
  <xsl:if test="not(local-name($type_node)='ReferenceType' or local-name($type_node)='PointerType')"> const&amp;</xsl:if>
  <xsl:value-of select="concat(' ', @name)"/>
 </xsl:for-each>
 <xsl:if test="$has_return_arg">
  <xsl:if test="not(@static='1') or not(empty(Argument)) or not($method/@static='1')">, </xsl:if>
  <xsl:choose>
   <xsl:when test="@returns=$const_ref_string_type_id or @returns=$ptr_const_char_type_id">
    <xsl:value-of select="'char*'"/> <!-- We need a non-const buffer to copy the result into -->
   </xsl:when>
   <xsl:otherwise>
    <xsl:call-template name="wrap_type"><xsl:with-param name="type_id" select="@returns"/></xsl:call-template>
    <xsl:value-of select="'&amp;'"/>
   </xsl:otherwise>
  </xsl:choose>
  <xsl:value-of select="' result'"/>
 </xsl:if>
 <!-- If any argument is a string, include a length argument -->
 <xsl:for-each select="Argument">
  <xsl:if test="@type=$const_ref_string_type_id or @type=$ptr_const_char_type_id">
   <xsl:value-of select="concat(', int ', @name, '_length')"/>
  </xsl:if>
 </xsl:for-each>
 <xsl:if test="$has_return_arg and (@returns=$const_ref_string_type_id or @returns=$ptr_const_char_type_id)">
  <xsl:value-of select="', int result_length'"/>
 </xsl:if>
 <xsl:value-of select="concat(') {', $newline, '    ')"/>
 <!-- Now the method body -->
 <xsl:choose>
  <xsl:when test="$has_return">
   <xsl:value-of select="'return '"/>
  </xsl:when>
  <xsl:when test="$has_return_arg and (@returns=$const_ref_string_type_id or @returns=$ptr_const_char_type_id)">
   <xsl:value-of select="'const char* result_chars = '"/>
  </xsl:when>
  <xsl:when test="$has_return_arg">
   <xsl:value-of select="'result = '"/>
  </xsl:when>
 </xsl:choose>
 <xsl:value-of select="concat('OpenMM_', $class_name, '_', @name, '(')"/>
 <xsl:if test="not(@static='1')">target</xsl:if>
 <xsl:for-each select="Argument">
  <xsl:if test="position() > 1 or not($method/@static='1')">, </xsl:if>
   <xsl:variable name="type_id" select="@type"/>
   <xsl:variable name="type_node" select="/GCC_XML/*[@id=$type_id]"/>
   <xsl:choose>
   <xsl:when test="@type=$const_ref_string_type_id or @type=$ptr_const_char_type_id">
    <xsl:value-of select="concat('makeString(', @name, ', ', @name, '_length).c_str()')"/>
   </xsl:when>
   <xsl:when test="local-name($type_node)='Enumeration'">
    <xsl:variable name="enum_class" select="/GCC_XML/*[@id=$type_node/@context]"/>
    <xsl:value-of select="concat('(OpenMM_', $enum_class/@name, '_', $type_node/@name, ') ', @name)"/>
   </xsl:when>
   <xsl:when test="local-name($type_node)='ReferenceType' and not(empty(/GCC_XML/Enumeration[@id=$type_node/@type]))">
    <xsl:variable name="enum_node" select="/GCC_XML/Enumeration[@id=$type_node/@type]"/>
    <xsl:variable name="enum_class" select="/GCC_XML/*[@id=$enum_node/@context]"/>
    <xsl:value-of select="concat('(OpenMM_', $enum_class/@name, '_', $enum_node/@name, '*) ', @name)"/>
   </xsl:when>
   <xsl:otherwise>
    <xsl:value-of select="@name"/>
   </xsl:otherwise>
  </xsl:choose>
 </xsl:for-each>
 <xsl:value-of select="');'"/>
 <xsl:if test="$has_return_arg and (@returns=$const_ref_string_type_id or @returns=$ptr_const_char_type_id)">
    copyAndPadString(result, result_chars, result_length);</xsl:if>
};
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
  <xsl:when test="$type_id=$vector_bond_type_id">
   <xsl:value-of select="'OpenMM_BondArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_2d_int_type_id">
   <xsl:value-of select="'OpenMM_2D_IntArray'"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_tortor_type_id">
   <xsl:value-of select="'OpenMM_3D_DoubleArray'"/>
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
   <xsl:value-of select="'int'"/>
  </xsl:when>
  <xsl:when test="$node/@context=$openmm_namespace_id">
   <xsl:value-of select="concat('OpenMM_', $node/@name)"/>
  </xsl:when>
  <xsl:otherwise>
   <xsl:value-of select="$node/@name"/>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<!-- Determine whether a type is one of the special handle types -->
<xsl:template name="is_handle_type">
 <xsl:param name="type_id"/>
 <xsl:variable name="node" select="/GCC_XML/*[@id=$type_id]"/>
 <xsl:choose>
  <xsl:when test="$type_id=$vec3_type_id"></xsl:when>
  <xsl:when test="$type_id=$vector_string_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_vec3_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_bond_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_2d_int_type_id">1</xsl:when>
  <xsl:when test="$type_id=$map_parameter_type_id">1</xsl:when>
  <xsl:when test="$type_id=$map_property_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_double_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_int_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_string_type_id">1</xsl:when>
  <xsl:when test="$type_id=$vector_tortor_type_id">1</xsl:when>
  <xsl:when test="local-name($node)='Class' and $node/@context=$openmm_namespace_id">1</xsl:when>
  <xsl:when test="local-name($node)='ReferenceType' or local-name($node)='PointerType'">
   <xsl:call-template name="is_handle_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='CvQualifiedType'">
   <xsl:call-template name="is_handle_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
  </xsl:when>
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
  <xsl:when test="local-name($node)='Enumeration' or (local-name($node)='ReferenceType' and not(empty(/GCC_XML/Enumeration[@id=$node/@type])))">
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
