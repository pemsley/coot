<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:strip-space elements="*"/>

<!-- replace all occurrences in a string-->
<xsl:template name="string-replace-all">
  <xsl:param name="string"/>
  <xsl:param name="replace"/>
  <xsl:param name="by"/>
  <xsl:choose>
    <xsl:when test="contains($string, $replace)">
      <xsl:value-of select="substring-before($string, $replace)"/>
      <xsl:value-of select="$by"/>
      <xsl:call-template name="string-replace-all">
        <xsl:with-param name="string" select="substring-after($string,$replace)"/>
        <xsl:with-param name="replace" select="$replace"/>
        <xsl:with-param name="by" select="$by"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$string"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- convert argument type to something easily readable -->
<xsl:template name="guess-argument-type">
  <xsl:param name="arg"/>
  <xsl:choose>
    <xsl:when test="contains($arg,'int')">
      <xsl:text>an integer number</xsl:text>
    </xsl:when>
    <xsl:when test="contains($arg,'long')">
      <xsl:text>an integer number</xsl:text>
    </xsl:when>
    <xsl:when test="contains($arg,'double')">
      <xsl:text>a number</xsl:text>
    </xsl:when>
    <xsl:when test="contains($arg,'float')">
      <xsl:text>a number</xsl:text>
    </xsl:when>
    <xsl:when test="contains($arg,'char')">
      <xsl:choose>
        <xsl:when test="contains($arg,'char *')">
          <xsl:text>a string</xsl:text>
        </xsl:when>
        <xsl:when test="$arg='char'">
          <xsl:text>a character</xsl:text>
        </xsl:when>
      </xsl:choose>
    </xsl:when>
    <xsl:otherwise> <!-- whatever without the " *" -->
      <xsl:text>a </xsl:text>
      <xsl:call-template name="string-replace-all">
        <xsl:with-param name="string" select="$arg"/>
        <xsl:with-param name="replace" select="' *'"/>
        <xsl:with-param name="by" select="''"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- doxy to texi-->
<xsl:output method="text"/>
<xsl:template match="/doxygen">

  <!-- generate menu-->
  <xsl:text>@menu&#xa;</xsl:text>
  <xsl:for-each select="child::*/child::sectiondef">
    <xsl:if test="@kind='user-defined' and
                  count(child::memberdef/child::briefdescription/child::para)&gt;0">
      <xsl:text>* </xsl:text>
      <xsl:value-of select="header"/><xsl:text>::&#xa;</xsl:text>
    </xsl:if>
    <xsl:if test="@kind='func' and
                  count(child::memberdef/child::briefdescription/child::para)&gt;0">
       <xsl:text>* Sectionless functions</xsl:text><xsl:text>::&#xa;</xsl:text>
    </xsl:if>
  </xsl:for-each>
  <xsl:text>@end menu&#xa;&#xa;</xsl:text>

  <!-- find sections -->
  <xsl:for-each select="child::*/child::sectiondef">
    <xsl:if test="(@kind='user-defined' or @kind='func') and
                  count(child::memberdef/child::briefdescription/child::para)&gt;0">
                      
      <xsl:if test="@kind='user-defined'">
        <xsl:text>@node </xsl:text>
        <xsl:value-of select="header"/><xsl:text>&#xa;</xsl:text>
        <xsl:text>@section </xsl:text>
        <xsl:value-of select="header"/><xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text>@menu</xsl:text>
        <xsl:text>&#xa;</xsl:text>
      </xsl:if>
      <xsl:if test="@kind='func'">
        <xsl:text>@node Sectionless functions&#xa;</xsl:text>
        <xsl:text>@section Sectionless functions&#xa;</xsl:text>
        <xsl:text>@menu</xsl:text>
        <xsl:text>&#xa;</xsl:text>
      </xsl:if>

      <!-- generate menu for functions -->
      <xsl:for-each select="child::memberdef">
        <xsl:if test="@kind='function' and
                      count(child::briefdescription/child::para)&gt;0">
          <xsl:variable name="FunctionName">
            <xsl:value-of select="name"/>
          </xsl:variable>

          <xsl:text>* </xsl:text>
          <xsl:call-template name="string-replace-all">
            <xsl:with-param name="string" select="$FunctionName"/>
            <xsl:with-param name="replace" select="'_'"/>
            <xsl:with-param name="by" select="'-'"/>
          </xsl:call-template>
          <xsl:text>::&#xa;</xsl:text>
        </xsl:if>
      </xsl:for-each>
      <xsl:text>@end menu</xsl:text>
      <xsl:text>&#xa;&#xa;</xsl:text>

      <!-- find functions-->
      <xsl:for-each select="child::memberdef">
      <xsl:if test="@kind='function' and
                   count(child::briefdescription/child::para)&gt;0">
        <xsl:variable name="FunctionName">
          <xsl:value-of select="name"/>
        </xsl:variable>

        <xsl:variable name="SchemeyName">
          <xsl:call-template name="string-replace-all">
            <xsl:with-param name="string" select="$FunctionName"/>
            <xsl:with-param name="replace" select="'_'"/>
            <xsl:with-param name="by" select="'-'"/>
          </xsl:call-template>
        </xsl:variable>

        <xsl:text>@node </xsl:text>
        <xsl:value-of select="$SchemeyName"/>
        <xsl:text>&#xa;</xsl:text>
        <xsl:text>@subsection </xsl:text>
        <xsl:value-of select="$SchemeyName"/>
        <xsl:text>&#xa;</xsl:text>
        <xsl:text>@deffn </xsl:text>
        <xsl:value-of select="@kind"/>
        <xsl:text> </xsl:text>
        <xsl:value-of select="$SchemeyName"/>


        <!-- parameters -->
        <!-- parameters for @deffn line-->
        <xsl:for-each select="child::param">
          <xsl:text> </xsl:text>
          <xsl:value-of select="declname"/>
        </xsl:for-each>

        <!-- parameters for listing -->
        <xsl:if test="count(child::param)>1">
        <xsl:text>&#xa;&#xa;</xsl:text>
          <xsl:text>Where: &#xa;@itemize @bullet&#xa;</xsl:text>
          <xsl:for-each select="child::param">
            <xsl:text>@item @emph{</xsl:text>
            <xsl:value-of select="declname"/>
            <xsl:text>} is </xsl:text>
            <xsl:variable name="Type">
              <xsl:value-of select="type"/>
            </xsl:variable>
            <xsl:call-template name="guess-argument-type">
              <xsl:with-param name="arg" select="$Type"/>
            </xsl:call-template>
            <xsl:text>&#xa;</xsl:text>
          </xsl:for-each>
          <xsl:text>@end itemize</xsl:text>
        </xsl:if>

        <xsl:if test="count(child::param)=1">
          <xsl:variable name="Type">
            <xsl:value-of select="child::param/child::type"/>
          </xsl:variable>
          <xsl:text>&#xa;&#xa;</xsl:text>
          <xsl:text>Where @emph{</xsl:text>
          <xsl:value-of select="child::param/child::declname"/>
            <xsl:text>} is </xsl:text>
            <xsl:call-template name="guess-argument-type">
              <xsl:with-param name="arg" select="$Type"/>
            </xsl:call-template>
        </xsl:if>

        <xsl:text>&#xa;&#xa;</xsl:text>

        <!-- descriptions -->
        <xsl:for-each select="child::briefdescription/child::para">
          <xsl:value-of select="."/>
          <xsl:text>&#xa;&#xa;</xsl:text>
        </xsl:for-each>

        <xsl:for-each select="detaileddescription/para/child::text()">
           <xsl:value-of select="."/>
           <xsl:text>&#xa;&#xa;</xsl:text>
        </xsl:for-each>

        <xsl:for-each select="detaileddescription/para/simplesect/para/child::text()">
           <xsl:text>Returns: </xsl:text>
           <xsl:value-of select="."/>
           <xsl:text>&#xa;&#xa;</xsl:text>
        </xsl:for-each>

        <xsl:for-each select="child::inbodydescription/child::para">
          <xsl:value-of select="."/>
          <xsl:text>&#xa;&#xa;</xsl:text>
        </xsl:for-each>

        <xsl:text>@end deffn&#xa;&#xa;</xsl:text>

      </xsl:if> <!-- widget function filtering -->
      </xsl:for-each>   <!--function-->

    </xsl:if>
  </xsl:for-each>   <!--section-->

</xsl:template>
</xsl:stylesheet>

